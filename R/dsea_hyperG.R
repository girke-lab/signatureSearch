##' Drug Set Enrichment Analysis (DSEA) with Hypergeometric Test
##' 
##' The \code{dsea_hyperG} function performs Drug Set Enrichment Analysis (DSEA)
##' based on the hypergeomtric distribution. In case of DSEA, the identifiers of
##' the top ranking drugs from a GESS result table are used. To use drug 
##' instead of gene labels for this test, the former are mapped to functional 
##' categories, including GO, KEGG or Mode of Action (MOA) categories, based on 
##' drug-target interaction annotations provided by databases such as DrugBank, 
##' ChEMBL, CLUE or STITCH. Currently, the MOA annotation used by this function 
##' are from the CLUE website (https://clue.io).
##' 
##' Compared to the related Target Set Enrichment Analysis (TSEA; see help
##' \code{tsea_dup_hyperG} or \code{tsea_mGSEA}), the DSEA approach has the
##' advantage that the drugs in the query test sets are usually unique allowing 
##' to use them without major modifications to the underlying statistical 
##' method(s).
##' 
##' The DSEA results stored in the \code{feaResult} object can be returned with 
##' the \code{result} method in tabular format, here \code{tibble}. The columns 
##' of this \code{tibble} are described in the help of the 
##' \code{\link{tsea_dup_hyperG}} function.
##' @param drugs character vector, query drug identifier set used for functional
##' enrichment testing. This can be the top ranking drugs from a GESS result.
##' @param type one of 'GO', 'KEGG' or 'MOA'
##' @param ont character(1). If type is `GO`, assign \code{ont} (ontology) one 
##' of `BP`,`MF`, `CC` or `ALL`. If type is 'KEGG', \code{ont} is ignored.
##' @param pvalueCutoff double, p-value cutoff to return only enrichment results
##' for drugs meeting a user definable confidence threshold
##' @param pAdjustMethod p-value adjustment method, 
#' one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'
##' @param qvalueCutoff double, qvalue cutoff, similar to \code{pvalueCutoff}
##' @param minGSSize integer, annotation categories with less than 
##' \code{minGSize} drugs annotated will be ignored by enrichment test. If type 
##' is 'MOA', it may be beneficial to set 'minGSSize' to lower values (e.g. 2) 
##' than for other functional annotation systems. This is because certain MOA 
##' categories contain only 2 drugs.
##' @param maxGSSize integer, annotation categories with more drugs annotated 
##' than \code{maxGSize} will be ignored by enrichment test.
##' @return \code{\link{feaResult}} object containing the enrichment results of 
##' functional categories (e.g. GO terms or KEGG pathways) ranked by the 
##' corresponding enrichment statistic.
##' @import org.Hs.eg.db
##' @importFrom magrittr %<>%
##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importMethodsFrom AnnotationDbi mget
##' @seealso \code{\link{feaResult}}, \code{\link{fea}},
##'          \code{\link[signatureSearchData]{GO_DATA_drug}}
##' @examples 
##' data(drugs10)
##' ## GO annotation system
##' # hyperG_res <- dsea_hyperG(drugs = drugs10, type = "GO", ont="MF")
##' # result(hyperG_res)
##' ## KEGG annotation system
##' #hyperG_k_res <- dsea_hyperG(drugs = drugs10, type = "KEGG", 
##' #                            pvalueCutoff = 1, qvalueCutoff = 1, 
##' #                            minGSSize = 10, maxGSSize = 500)
##' #result(hyperG_k_res) 
##' @export
dsea_hyperG <- function(drugs,
                        type="GO",
                        ont="BP",
                        pvalueCutoff=0.05,
                        pAdjustMethod="BH",
                        qvalueCutoff=0.2,
                        minGSSize=10,
                        maxGSSize=500) {
  drugs <- as.character(unique(tolower(drugs)))
  if(type == "GO"){
    ont %<>% toupper
    ont <- match.arg(ont, c("BP", "CC", "MF", "ALL"))
    
    # GO_DATA_drug <- get_GO_data_drug(OrgDb = "org.Hs.eg.db", 
    #                                  ont, keytype="SYMBOL")
    # download GO_DATA_drug.rds from AnnotationHub to save time by avoiding 
    # builing it from scratch
    GO_DATA_drug <- suppressMessages(ah[["AH69087"]])
    
    # get all the drugs in the corresponding annotation system as universe
    ext2path <- get("EXTID2PATHID", envir = GO_DATA_drug)
    universe = names(ext2path)
    
    res <- enricher_internal(gene=drugs,
                             pvalueCutoff=pvalueCutoff,
                             pAdjustMethod=pAdjustMethod,
                             universe = universe,
                             qvalueCutoff = qvalueCutoff,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             USER_DATA = GO_DATA_drug)
    
    if (is.null(res))
      return(res)
    
    og(res) <- get_organism(OrgDb="org.Hs.eg.db")
    ont(res) <- ont
    tg(res) <- NULL
    if (ont == "ALL") {
      res <- add_GO_Ontology(res, GO_DATA_drug)
    }
    return(res)
  }
  if(type == "KEGG"){
    species <- organismMapper("hsa")
    KEGG_DATA_drug <- prepare_KEGG_drug(species, "KEGG", keyType="kegg")
    # get all the drugs in the corresponding annotation system as universe
    ext2path <- get("EXTID2PATHID", envir = KEGG_DATA_drug)
    universe = names(ext2path)
    res <- enricher_internal(drugs,
                             pvalueCutoff  = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             universe      = universe,
                             minGSSize     = minGSSize,
                             maxGSSize     = maxGSSize,
                             qvalueCutoff  = qvalueCutoff,
                             USER_DATA = KEGG_DATA_drug)
    if (is.null(res))
      return(res)
    tg(res) <- NULL
    ont(res) <- "KEGG"
    og(res) <- species
    return(res)
  }
  if(type == "MOA"){
      data("clue_moa_list", envir = environment())
      moa_list <- clue_moa_list
      MOA_DATA <- get_MOA_data(moa_list, keytype="drug_name")
      # get all the drugs in the corresponding annotation system as universe
      ext2path <- get("EXTID2PATHID", envir = MOA_DATA)
      universe = names(ext2path)
      
      res <- enricher_internal(gene=drugs,
                               pvalueCutoff=pvalueCutoff,
                               pAdjustMethod=pAdjustMethod,
                               universe = universe,
                               qvalueCutoff = qvalueCutoff,
                               minGSSize = minGSSize,
                               maxGSSize = maxGSSize,
                               USER_DATA = MOA_DATA)
      
      if (is.null(res))
          return(res)
      
      og(res) <- "Homo sapiens"
      ont(res) <- "MOA"
      tg(res) <- NULL
      return(res)
  }
}
