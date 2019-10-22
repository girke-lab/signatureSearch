#' The \code{dsea_GSEA} function performs Drug Set Enrichment Analysis (DSEA)
#' with the GSEA algorithm from Subramanian et al. (2005). In case of DSEA, drug
#' identifiers combined with their ranking scores of an upstream GESS method are
#' used, such as the NCS values from the LINCS method. To use drug instead of
#' gene labels for GSEA, the former are mapped to functional categories, 
#' including GO or KEGG, based on drug-target interaction
#' annotations provided by databases such as DrugBank, ChEMBL, CLUE or
#' STITCH. 
#' 
#' The DSEA results stored in the \code{feaResult} object can be returned with 
#' the \code{result} method in tabular format, here \code{tibble}. The columns 
#' of this \code{tibble} are described in the help of the 
#' \code{\link{tsea_mGSEA}} function.
#' 
#' @title Drug Set Enrichment Analysis (DSEA) with GSEA Algorithm
#' @param drugList named numeric vector, where the names represent drug labels 
#' and the numeric component scores. This can be all drugs of a GESS result that
#' are ranked by GESS scores, such as NCSs of the LINCS method. Note, drugs
#' with scores of zero are ignored by this method.
#' 
#' @param type one of 'GO', 'KEGG' or 'MOA'
#' @param ont character(1). If type is `GO`, assign \code{ont} (ontology) one 
#' of `BP`,`MF`, `CC` or `ALL`. If type is 'KEGG', \code{ont} is ignored.
#' @param exponent integer value used as exponent in GSEA algorithm. It defines
#' the weight of the items in the item set \emph{S}. Note, in DSEA the items 
#' are drug labels, while it is gene labels in the original GSEA.
#' @param nPerm integer defining the number of permutation iterations for 
#' calculating p-values
#' @param minGSSize integer, annotation categories with less than 
#' \code{minGSize} drugs annotated will be ignored by enrichment test. If type 
#' is 'MOA', it may be beneficial to set 'minGSSize' to lower values (e.g. 2) 
#' than for other functional annotation systems. This is because certain MOA 
#' categories contain only 2 drugs.
#' @param maxGSSize integer, annotation categories with more drugs annotated 
#' than \code{maxGSize} will be ignored by enrichment test.
#' @param pvalueCutoff double, p-value cutoff to return only enrichment results
#' for drugs meeting a user definable confidence threshold
#' @param pAdjustMethod p-value adjustment method, 
#' one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'
#' @return \code{\link{feaResult}} object containing the enrichment results of 
#' functional categories (e.g. GO terms or KEGG pathways) ranked by the 
#' corresponding enrichment statistic.
#' @seealso \code{\link{feaResult}}, \code{\link{fea}},
#'          \code{\link[signatureSearchData]{GO_DATA_drug}}
#' @references 
#' GSEA algorithm: 
#' Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., 
#' Gillette, M. A., Mesirov, J. P. (2005). Gene set enrichment analysis: a 
#' knowledge-based approach for interpreting genome-wide expression profiles. 
#' Proceedings of the National Academy of Sciences of the United States of
#' America, 102(43), 15545-15550. URL: https://doi.org/10.1073/pnas.0506580102
#' @examples 
#' data(drugs10)
#' dl <- c(rev(seq(0.1, 0.5, by=0.05)), 0)
#' names(dl)=drugs10
#' ## KEGG annotation system
#' #gsea_k_res <- dsea_GSEA(drugList=dl, type="KEGG", exponent=1, nPerm=100, 
#' #                        pvalueCutoff=0.5, minGSSize=2)
#' #result(gsea_k_res)
#' @export
dsea_GSEA <- function(drugList, 
                      type = "GO", 
                      ont           = "BP", 
                      exponent      = 1,
                      nPerm         = 1000,
                      minGSSize     = 10,
                      maxGSSize     = 500,
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH") {
  
  ont %<>% toupper
  ont <- match.arg(ont, c("BP", "CC", "MF", "ALL"))
  names(drugList) <- tolower(names(drugList))
  drugList <- drugList[drugList>0]
  # check whether there are duplicated names in drugList
  if(sum(duplicated(names(drugList))) > 0)
    stop("The names of the query drug list for dsea_GSEA need to be unique!")
  
  if(type=="GO"){
    # GO_DATA_drug <- get_GO_data_drug(OrgDb = "org.Hs.eg.db", 
    #                                  ont, keytype="SYMBOL")
    # download GO_DATA_drug.rds 
    eh <- suppressMessages(ExperimentHub())
    GO_DATA_drug <- suppressMessages(eh[["EH3232"]])
    
    res <-  GSEA_internal(geneList = drugList,
                          exponent = exponent,
                          nPerm = nPerm,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          USER_DATA = GO_DATA_drug,
                          seed = FALSE)
    
    if (is.null(res))
      return(res)
    drugs(res) <- names(drugList)
    tg(res) <- NULL
    og(res) <- get_organism(OrgDb = "org.Hs.eg.db")
    ont(res) <- ont
    
    if (ont == "ALL") {
      res <- add_GO_Ontology(res, GO_DATA_drug)
    } 
    return(res)
  }
  
  if(type=="KEGG"){
    species <- organismMapper("hsa")
    KEGG_DATA_drug <- prepare_KEGG_drug(species, "KEGG", keyType="kegg")
    
    res <-  GSEA_internal(geneList = drugList,
                          exponent = exponent,
                          nPerm = nPerm,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          USER_DATA = KEGG_DATA_drug,
                          seed = FALSE)
    
    if (is.null(res))
      return(res)
    drugs(res) <- names(drugList)
    tg(res) <- NULL
    og(res) <- species
    ont(res) <- "KEGG"
    return(res)
  }
  if(type == "MOA"){
      data("clue_moa_list", envir = environment())
      moa_list <- clue_moa_list
      MOA_DATA <- get_MOA_data(moa_list, keytype="drug_name")
      res <- GSEA_internal(geneList = drugList,
                           exponent = exponent,
                           nPerm = nPerm,
                           minGSSize = minGSSize,
                           maxGSSize = maxGSSize,
                           pvalueCutoff = pvalueCutoff,
                           pAdjustMethod = pAdjustMethod,
                           USER_DATA = MOA_DATA,
                           seed = FALSE)
      
      if (is.null(res))
          return(res)
      drugs(res) <- names(drugList)
      tg(res) <- NULL
      og(res) <- "Homo sapiens"
      ont(res) <- "MOA"
      return(res)
  }
}