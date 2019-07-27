#' The classical GSEA method is used to do enrichment analysis on a scored
#' ranked drug list after mapping drugs to functional categories via 
#' drug-target links in DrugBank, CLUE and STITCH databases for GO and KEGG 
#' pathway enrichment. It can also be used to get the enriched MOAs of a 
#' scored ranked drug list. 
#' 
#' In the GESS/FEA workflow, the ranked drug list consists of all drugs in the 
#' GESS result. The similarity scores of the corresponding GESS method can be 
#' used for ranking the drugs as it is required by the GSEA algorithm. The drugs
#' with zero scores are excluded. 
#'
#' The description of columns in the result table can be found at 
#' \code{\link{tsea_mGSEA}} function.
#' 
#' @title GSEA method for DSEA
#' @param drugList scored ranked list of all drugs in the GESS result. 
#' The similarity scores of the corresponding GESS method can be used for 
#' ranking the drugs as it is required by the GSEA algorithm. 
#' The drugs with zero scores are excluded.
#' @param type one of "GO", "KEGG" or "MOA"
#' @param ont one of "BP", "MF", "CC" or "ALL" if type is "GO"
#' @param exponent weight of each step
#' @param nPerm permutation numbers
#' @param minGSSize minimal size of drug sets annotated by ontology term 
#' after drug to functional category mappings. If type is "MOA", it is
#' recommended to set minGSSize as 2 since some MOA categories only contain 2 
#' drugs
#' @param maxGSSize maximal size of drug sets annotated for testing
#' @param pvalueCutoff pvalue Cutoff
#' @param pAdjustMethod p value adjustment method,
#' one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @return \code{\link{feaResult}} object, 
#' represents enriched functional categories.
#' @seealso \code{\link{feaResult}}, \code{\link{fea}},
#'          \code{\link[signatureSearchData]{dtlink_db_clue_sti}}
#' @references 
#' GSEA algorithm: 
#' Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., 
#' Gillette, M. A., … Mesirov, J. P. (2005). Gene set enrichment analysis: a 
#' knowledge-based approach for interpreting genome-wide expression profiles. 
#' Proceedings of the National Academy of Sciences of the United States of
#' America, 102(43), 15545–15550. \url{https://doi.org/10.1073/pnas.0506580102}
#' @examples 
#' db_path <- system.file("extdata", "sample_db.h5", 
#'                        package = "signatureSearch")
#' library(signatureSearchData)
#' sample_db <- readHDF5chunk(db_path, colindex=1:100)
#' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
#' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
#' query = as.numeric(query_mat); names(query) = rownames(query_mat)
#' upset <- head(names(query[order(-query)]), 150)
#' downset <- tail(names(query[order(-query)]), 150)
#' qsig_lincs <- qSig(query = list(upset=upset, downset=downset), 
#'                    gess_method = "LINCS", refdb = db_path)
#' lincs <- gess_lincs(qsig_lincs, sortby="NCS", tau=FALSE)
#' dl <- abs(result(lincs)$NCS); names(dl) <- result(lincs)$pert
#' dl <- dl[dl>0]
#' dl <- dl[!duplicated(names(dl))]
#' ## GO annotation system
#' \dontrun{
#' gsea_res <- dsea_GSEA(drugList=dl, type="GO", ont="MF", exponent=1, 
#'                       nPerm=1000, pvalueCutoff=0.2, minGSSize=5)
#'                       result(gsea_res)
#' }
#' ## KEGG annotation system
#' gsea_k_res <- dsea_GSEA(drugList=dl, type="KEGG", exponent=1, nPerm=1000, 
#'                         pvalueCutoff=0.5, minGSSize=5)
#' result(gsea_k_res)
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
  # check whether there are duplicated names in drugList
  if(sum(duplicated(names(drugList))) > 0)
    stop("The names of the query drug list for dsea_GSEA need to be unique!")
  
  if(type=="GO"){
    # GO_DATA_drug <- get_GO_data_drug(OrgDb = "org.Hs.eg.db", 
    #                                  ont, keytype="SYMBOL")
    # download GO_DATA_drug.rds 
    GO_DATA_drug <- suppressMessages(ah[["AH69087"]])
    
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
    res@drugs <- names(drugList)
    res@targets <- NULL
    res@organism <- get_organism(OrgDb = "org.Hs.eg.db")
    res@ontology <- ont
    
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
    res@drugs <- names(drugList)
    res@targets <- NULL
    res@organism <- species
    res@ontology <- "KEGG"
    return(res)
  }
  if(type == "MOA"){
      moa_list <- readRDS(system.file("extdata", "clue_moa_list.rds", 
                                      package="signatureSearch"))
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
      res@drugs <- names(drugList)
      res@targets <- NULL
      res@organism <- "Homo sapiens"
      res@ontology <- "MOA"
      return(res)
  }
}



