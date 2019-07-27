##' hyperG method for DSEA
##' 
##' The hypergeometric test is used to do enrichment analysis on a drug set 
##' after mapping drugs to functional categories via drug-target links in
##' DrugBank, CLUE and STITCH databases for GO and KEGG pathway enrichment. 
##' It can also be used to get the enriched MOAs of a drug set. 
##' The drugs MOA annotation come from CLUE website.
##' 
##' The drug sets can be directly used for GO/KEGG pathways enrichment testing 
##' by changing the mappings in the reference database from target-to-functional
##' category mappings to drug-to-functional category mappings. 
##' The latter can be generated based on the drug-target information provided 
##' by DrugBank or related databases. As a result, one can perform the FEA on 
##' ranked drug lists directly. This DSEA approach has the advantage that the 
##' drugs in the query test sets are usually unique allowing to use them 
##' without any changes for functional enrichment methods. 
##' 
##' So the hypergeometric test is directly used when the query is a vector of 
##' drugs, e.g., top ranking drugs from GESS result since the functional
##' categories also contains drug sets after gene to drug mappings.
##'
##' The description of columns in the result table can be found at 
##' \code{\link{tsea_dup_hyperG}} function.
##' 
##' @param drugs query drug set used to do DSEA. 
##' Can be top ranking drugs in GESS result. 
##' @param type one of "GO", "KEGG" or "MOA"
##' @param ont One of "MF", "BP", "CC" or "ALL" if type is "GO"
##' @param pvalueCutoff Cutoff value of p value.
##' @param pAdjustMethod p value adjustment method,
##' one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param qvalueCutoff q value cutoff
##' @param minGSSize minimal size of drug sets annotated by ontology term 
##' after drug to functional category mappings. If type is "MOA", it is
##' recommended to set minGSSize as 2 since some MOA categories only contain 2 
##' drugs.
##' @param maxGSSize maximal size of drug sets annotated for testing
##' @return \code{\link{feaResult}} object, 
##' represents enriched functional categories.
##' @import org.Hs.eg.db
##' @importFrom magrittr %<>%
##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importMethodsFrom AnnotationDbi mget
##' @seealso \code{\link{feaResult}}, \code{\link{fea}},
##'          \code{\link[signatureSearchData]{dtlink_db_clue_sti}}
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
