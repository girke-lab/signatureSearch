#' @rdname fea
#' @description 
#' The TSEA with meanAbs (\code{tsea_mabs}) method is a simple but effective 
#' functional enrichment statistic (Fang et al., 2012). As required for TSEA, 
#' it supports query label sets (here for target proteins/genes) with 
#' duplications by transforming them to score ranked label lists and then 
#' calculating mean absolute scores of 
#' labels in label set \eqn{S}.
#' @details 
#' The input for the mabs method is \eqn{L_tar}, the same as for mGSEA. In this
#' enrichment statistic, \eqn{mabs(S)}, of a label (e.g. gene/protein) set
#' \eqn{S} is calculated as mean absolute scores of the labels in \eqn{S}. In
#' order to adjust for size variations in label set \eqn{S}, 1000 random
#' permutations of \eqn{L_tar} are performed to determine \eqn{mabs(S,pi)}.
#'Subsequently, \eqn{mabs(S)} is normalized by subtracting the median of the
#' \eqn{mabs(S,pi)} and then dividing by the standard deviation of
#' \eqn{mabs(S,pi)} yielding the normalized scores \eqn{Nmabs(S)}. Finally, the
#' portion of \eqn{mabs(S,pi)} that is greater than \eqn{mabs(S)} is used as
#' nominal p-value (Fang et al., 2012). The resulting nominal p-values are
#' adjusted for multiple hypothesis testing using the Benjamini-Hochberg method.
#' @examples 
#' 
#' ############# MeanAbs method ##############
#' ## GO annotation system
#' # res1 <- tsea_mabs(drugs=drugs10, type="GO", ont="MF", nPerm=1000, 
#' #                   pvalueCutoff=0.05, minGSSize=5)
#' # result(res1)
#' ## KEGG annotation system
#' # res2 <- tsea_mabs(drugs=drugs10, type="KEGG", nPerm=1000, 
#' #                   pvalueCutoff=0.05, minGSSize=5)
#' # result(res2)
#' ## Reactome annotation system
#' # res3 <- tsea_mabs(drugs=drugs10, type="Reactome", pvalueCutoff=1)
#' # result(res3)
#' @export
#' 
tsea_mabs <- function(drugs, 
                      type="GO", ont="MF", 
                      nPerm=1000,  
                      pAdjustMethod="BH", pvalueCutoff=0.05,
                      minGSSize=5, maxGSSize=500, 
                      dt_anno="all", readable=FALSE){
  if(!any(type %in% c("GO", "KEGG", "Reactome"))){
    stop('"type" argument needs to be one of "GO", "KEGG" or "Reactome"')
  }
  drugs <- unique(tolower(drugs))
  targets <- get_targets(drugs, database = dt_anno)
  gnset <- na.omit(unlist(lapply(targets$t_gn_sym, function(i) 
    unlist(strsplit(as.character(i), split = "; ")))))
  # give scores to gnset
  tar_tab <- table(gnset)
  tar_dup <- as.numeric(tar_tab); names(tar_dup) <- names(tar_tab)
  tar_weight <- sort(tar_dup/sum(tar_dup), decreasing = TRUE)
  
  if(type=="GO"){
    # Get universe genes in GO annotation system
    universe <- univ_go
    tar_diff <- setdiff(universe, gnset)
    tar_diff_weight <- rep(0, length(tar_diff))
    names(tar_diff_weight) <- tar_diff
    tar_total_weight <- c(tar_weight, tar_diff_weight)
    mabsgo <- mabsGO(geneList = tar_total_weight, OrgDb = "org.Hs.eg.db", 
                     ont = ont, keyType = "SYMBOL", nPerm = nPerm, 
                     minGSSize = minGSSize, maxGSSize = maxGSSize, 
                     pvalueCutoff = pvalueCutoff, pAdjustMethod=pAdjustMethod)
    if(is.null(mabsgo))
        return(NULL)
    drugs(mabsgo) <- drugs
    return(mabsgo)
  }

  # convert gnset symbol to entrez
  OrgDb <- load_OrgDb("org.Hs.eg.db")
  gnset_map <- suppressMessages(AnnotationDbi::select(OrgDb, keys = gnset, 
                                  keytype = "SYMBOL", columns = "ENTREZID"))
  gnset_entrez <- as.character(na.omit(gnset_map$ENTREZID))
  # give scores to gnset_entrez
  tar_tab <- table(gnset_entrez)
  tar_dup <- as.numeric(tar_tab); names(tar_dup) <- names(tar_tab)
  tar_weight <- sort(tar_dup/sum(tar_dup), decreasing = TRUE)
  
  if(type=="KEGG"){
    # Get universe genes in KEGG annotation system
    KEGG_DATA <- prepare_KEGG(species="hsa", "KEGG", keyType="kegg")
    keggterms <- get("PATHID2EXTID", KEGG_DATA)
    universe <- unique(unlist(keggterms))
    
    tar_diff <- setdiff(universe, gnset_entrez)
    tar_diff_weight <- rep(0, length(tar_diff))
    names(tar_diff_weight) <- tar_diff
    tar_total_weight <- c(tar_weight, tar_diff_weight)

    mabs_res <- mabsKEGG(geneList=tar_total_weight, organism='hsa', 
                       keyType='kegg', nPerm = nPerm, 
                       minGSSize = minGSSize, maxGSSize=maxGSSize, 
                       pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod,
                       readable=readable)
  }
  
  if(type=="Reactome"){
    # Get universe genes in Reactome annotation system
    Reactome_DATA <- get_Reactome_DATA(organism="human")
    raterms <- get("PATHID2EXTID", Reactome_DATA)
    universe <- unique(unlist(raterms))
    
    tar_diff <- setdiff(universe, gnset_entrez)
    tar_diff_weight <- rep(0, length(tar_diff))
    names(tar_diff_weight) <- tar_diff
    tar_total_weight <- c(tar_weight, tar_diff_weight)
    
    mabs_res <- mabsReactome(geneList=tar_total_weight, organism='human', 
                             nPerm = nPerm, 
                             minGSSize = minGSSize, maxGSSize=maxGSSize, 
                             pvalueCutoff=pvalueCutoff, pAdjustMethod=pAdjustMethod,
                             readable=readable)
  }
  
  if(is.null(mabs_res))
    return(NULL)
  drugs(mabs_res) <- drugs
  return(mabs_res)
}