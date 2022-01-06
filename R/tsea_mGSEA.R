#' @rdname fea
#' @description 
#' The TSEA with mGSEA algorithm (\code{tsea_mGSEA} function) performs a 
#' Modified Gene Set Enrichment Analysis (mGSEA) that supports test sets 
#' (e.g. genes or protein IDs) with duplications. The duplication support is 
#' achieved by a weighting method for duplicated items, where the weighting is 
#' proportional to the frequency of the items in the test set.
#' @details 
#' The original GSEA method proposed by Subramanian et at., 2005 uses 
#' predefined gene sets \eqn{S} defined by functional annotation systems 
#' such as GO and KEGG. The goal is to determine whether the genes in \eqn{S}
#' are randomly distributed throughout a ranked test gene list \eqn{L} 
#' (e.g. all genes ranked by log2 fold changes) or enriched at the top or 
#' bottom of the test list. This is expressed by an 
#' Enrichment Score (\eqn{ES}) reflecting the degree to which a set \eqn{S} 
#' is overrepresented at the extremes of \eqn{L}. 
#' 
#' For TSEA, the query is a target protein set where duplicated entries need to 
#' be maintained. To perform GSEA with duplication support, here referred to as 
#' mGSEA, the target set is transformed to a score ranked target list 
#' \eqn{L_tar} of all targets provided by the 
#' corresponding annotation system. For each target in the query target set, 
#' its frequency is divided by the number of targets in the target set, 
#' which is the weight of that target. 
#' For targets present in the annotation system but absent in the 
#' target set, their scores are set to 0. Thus, every target in the annotation 
#' system will be assigned a score and then sorted decreasingly to obtain
#' \eqn{L_tar}.
#' 
#' In case of TSEA, the original GSEA method cannot be used directly since a 
#' large portion of zeros exists in \eqn{L_tar}. If the scores of the genes in 
#' set \eqn{S} are all zeros, \eqn{N_R} (sum of scores of genes in set 
#' \eqn{S}) will be zero, which cannot be used as the denominator. 
#' In this case, \eqn{ES} is set to -1. If only some genes in set \eqn{S}
#' have scores of zeros then \eqn{N_R} is set to a larger number to decrease 
#' the weight of the genes in \eqn{S} that have non-zero scores.
#'  
#' The reason for this modification is that if only one gene in gene set 
#' \eqn{S} has a non-zero score and this gene ranks high in \eqn{L_tar}, 
#' the weight of this gene will be 1 resulting in an \eqn{ES(S)} close to 1. 
#' Thus, the original GSEA method will score the gene set \eqn{S} as 
#' significantly enriched. However, this is undesirable because in this 
#' example only one gene is shared among the target set and the gene set 
#' \eqn{S}. Therefore, giving small weights (lowest non-zero score in \eqn{L_tar}) 
#' to genes in \eqn{S} that have zero scores could decrease the weight of the
#' genes in \eqn{S} that have non-zero scores, thereby decreasing the false 
#' positive rate. To favor truly enriched functional categories (gene set \eqn{S}) 
#' at the top of \eqn{L_tar}, only gene sets with positive \eqn{ES} are selected.
#' 
#' @param nPerm integer defining the number of permutation iterations for 
#' calculating p-values
#' @param exponent integer value used as exponent in GSEA algorithm. It defines
#' the weight of the items in the item set \eqn{S}.
#' @param verbose TRUE or FALSE, print message or not
#' @examples 
#' 
#' ############# TSEA mGSEA method ############
#' ## GO annotation system
#' # res1 <- tsea_mGSEA(drugs=drugs10, type="GO", ont="MF", exponent=1, 
#' #                    nPerm=1000, pvalueCutoff=1, minGSSize=5)
#' # result(res1)
#  ## KEGG annotation system
#' # res2 <- tsea_mGSEA(drugs=drugs10, type="KEGG", exponent=1, 
#' #                    nPerm=100, pvalueCutoff=1, minGSSize=5)
#' # result(res2)
#' ## Reactome annotation system
#' # res3 <- tsea_mGSEA(drugs=drugs10, type="Reactome", pvalueCutoff=1)
#' # result(res3)
#' @export
tsea_mGSEA <- function(drugs, type="GO", ont="MF", nPerm=1000, exponent=1, 
                       pAdjustMethod="BH", pvalueCutoff=0.05,
                       minGSSize=5, maxGSSize=500, verbose=FALSE,
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
    gsego <- gseGO2(geneList = tar_total_weight, OrgDb = "org.Hs.eg.db", 
                    ont = ont, keyType = "SYMBOL", nPerm = nPerm, 
                    minGSSize = minGSSize, maxGSSize = maxGSSize, 
                    exponent = exponent, nproc=1, verbose=verbose, 
                    pvalueCutoff = pvalueCutoff, pAdjustMethod=pAdjustMethod)
    if(is.null(gsego))
        return(NULL)
    drugs(gsego) <- drugs
    return(gsego)
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
    
    gse_res <- gseKEGG2(geneList=tar_total_weight, organism='hsa', keyType='kegg',
                        nPerm=nPerm, exponent=exponent,
                        minGSSize = minGSSize, maxGSSize=maxGSSize, 
                        pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod,
                        verbose = verbose, readable=readable)
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
    
    gse_res <- gseReactome(geneList=tar_total_weight, organism='human', 
                           exponent=exponent, nPerm=nPerm,
                           minGSSize = minGSSize, maxGSSize=maxGSSize, 
                           pvalueCutoff=pvalueCutoff, pAdjustMethod= pAdjustMethod,
                           verbose = verbose, readable=readable)
  }
  
  if(!is.null(gse_res))
    drugs(gse_res) <- drugs
  return(gse_res)
}
