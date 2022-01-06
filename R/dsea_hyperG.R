##' @rdname fea
##' @description 
##' The DSEA with Hypergeometric Test (\code{dsea_hyperG}) performs DSEA
##' based on the hypergeometric distribution. In case of DSEA, the identifiers of
##' the top ranking drugs from a GESS result table are used. To use drug 
##' instead of gene labels for this test, the former are mapped to functional 
##' categories, including GO, KEGG or Mode of Action (MOA) categories, based on 
##' drug-target interaction annotations provided by databases such as DrugBank, 
##' ChEMBL, CLUE or STITCH. Currently, the MOA annotation used by this function 
##' are from the CLUE website (https://clue.io).
##' 
##' Compared to the related Target Set Enrichment Analysis (TSEA), the DSEA 
##' approach has the advantage that the drugs in the query test sets are usually
##' unique allowing to use them without major modifications to the 
##' underlying statistical method(s).
##' @importFrom magrittr %<>%
##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importMethodsFrom AnnotationDbi mget
##' @examples 
##' 
##' ############### DSEA Hypergeometric Test ###########
##' ## GO annotation system
##' # hyperG_res <- dsea_hyperG(drugs=drugs10, type="GO", ont="MF")
##' # result(hyperG_res)
##' ## KEGG annotation system
##' # hyperG_k_res <- dsea_hyperG(drugs=drugs10, type="KEGG", 
##' #                             pvalueCutoff=1, qvalueCutoff=1, 
##' #                             minGSSize=10, maxGSSize=500)
##' # result(hyperG_k_res) 
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
    
    # get all the drugs in the corresponding annotation system as universe
    ext2path <- get("EXTID2PATHID", envir = GO_DATA_drug)
    universe <- names(ext2path)
    
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
    res <- select_ont(res, ont, GO_DATA_drug)
    og(res) <- get_organism(OrgDb = "org.Hs.eg.db")
    ont(res) <- ont
    tg(res) <- NULL
    return(res)
  }
  if(type == "KEGG"){
    species <- organismMapper("hsa")
    KEGG_DATA_drug <- prepare_KEGG_drug(species, "KEGG", keyType="kegg")
    # get all the drugs in the corresponding annotation system as universe
    ext2path <- get("EXTID2PATHID", envir = KEGG_DATA_drug)
    universe <- names(ext2path)
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
      universe <- names(ext2path)
      
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
