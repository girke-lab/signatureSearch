#' @rdname fea
#' @description 
#' The Drug Set Enrichment Analysis (DSEA) with GSEA algorithm 
#' (\code{dsea_GSEA} function) performs DSEA
#' with the GSEA algorithm from Subramanian et al. (2005). In case of DSEA, drug
#' identifiers combined with their ranking scores of an upstream GESS method are
#' used, such as the NCS values from the LINCS method. To use drug instead of
#' gene labels for GSEA, the former are mapped to functional categories, 
#' including GO or KEGG, based on drug-target interaction annotations provided 
#' by databases such as DrugBank, ChEMBL, CLUE or STITCH.
#' @examples 
#' data(drugs10)
#' 
#' ############ DSEA GSEA method ############
#' dl <- c(rev(seq(0.1, 0.5, by=0.05)), 0)
#' names(dl)=drugs10
#' ## KEGG annotation system
#' # gsea_k_res <- dsea_GSEA(drugList=dl, type="KEGG", exponent=1, nPerm=100, 
#' #                         pvalueCutoff=0.5, minGSSize=2)
#' # result(gsea_k_res)
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
    res <- select_ont(res, ont, GO_DATA_drug)
    drugs(res) <- names(drugList)
    tg(res) <- NULL
    og(res) <- get_organism(OrgDb = "org.Hs.eg.db")
    ont(res) <- ont
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
