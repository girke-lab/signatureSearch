##' hyperG Enrichment Method
##' 
##' This function uses hypergeometric test to do enrichment analysis on a 
##' drug set after mapping drugs to functional categories via drug-target 
##' links in DrugBank, CLUE and STITCH databases for GO and KEGG annotation 
##' system. It can also be used to get the enriched MOAs of a drug set if 
##' type is set as 'MOA'. The drugs MOA annotation is from CLUE website 
##' (\url{https://clue.io/}).
##' 
##' The drug set can be directly used to perform enrichment testing 
##' by changing the mappings in the annotation system from gene-to-functional
##' category mappings to drug-to-functional category mappings. 
##' The latter can be generated based on the drug-target information provided 
##' by DrugBank or related databases. As a result, one can perform the FEA on 
##' ranked drug list directly. This DSEA approach has the advantage that the 
##' drugs in the query test sets are usually unique allowing to use them 
##' without any changes for functional enrichment methods. 
##' 
##' Here, the hypergeometric test is directly used on a query drug test set,
##' e.g., top ranking drugs in GESS result since the functional
##' annotation system also contains drug sets after gene to drug mappings.
##' 
##' Note, description of the columns in the result table can be found at 
##' \code{\link{tsea_dup_hyperG}} function.
##' 
##' @param drugs character vector, query drug set used for functional 
##' enrichment. Can be top ranking drugs in the GESS result. 
##' @param type one of 'GO', 'KEGG' or 'MOA'
##' @param ont character(1). If type is `GO`, set ontology as one of `BP`,`MF`,
#' `CC` or `ALL`. If type is 'KEGG', it is ignored.
##' @param pvalueCutoff double, p-value cutoff
##' @param pAdjustMethod p-value adjustment method, 
#' one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'
##' @param qvalueCutoff double, qvalue cutoff
##' @param minGSSize integer, minimum size of each drug set in annotation 
##' system after drug to functional category mappings. If type is 'MOA', it is
##' recommended to set 'minGSSize' as 2 since some MOA categories only contain 
##' 2 drugs.
##' @param maxGSSize integer, maximum size of each drug set in annotation 
##' system 
##' @return \code{\link{feaResult}} object, the result table contains the
#' enriched functional categories (e.g. GO terms or KEGG pathways) ranked by 
#' the corresponding enrichment statistic.
##' @import org.Hs.eg.db
##' @importFrom magrittr %<>%
##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importMethodsFrom AnnotationDbi mget
##' @seealso \code{\link{feaResult}}, \code{\link{fea}},
##'          \code{\link[signatureSearchData]{goAnno_drug}}
##' @examples 
##' data(drugs)
##' ## GO annotation system
##' # hyperG_res <- dsea_hyperG(drugs = drugs, type = "GO", ont="MF")
##' # result(hyperG_res)
##' ## KEGG annotation system
##' hyperG_k_res <- dsea_hyperG(drugs = drugs, type = "KEGG", 
##'                             pvalueCutoff = 1, qvalueCutoff = 1, 
##'                             minGSSize = 10, maxGSSize = 2000)
##' result(hyperG_k_res) 
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
    
    res@organism <- get_organism(OrgDb="org.Hs.eg.db")
    res@ontology <- ont
    res@targets <- NULL
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
    res@targets <- NULL
    res@ontology <- "KEGG"
    res@organism <- species
    return(res)
  }
  if(type == "MOA"){
      moa_list <- readRDS(system.file("extdata", "clue_moa_list.rds", 
                                      package="signatureSearch"))
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
      
      res@organism <- "Homo sapiens"
      res@ontology <- "MOA"
      res@targets <- NULL
      return(res)
  }
}
