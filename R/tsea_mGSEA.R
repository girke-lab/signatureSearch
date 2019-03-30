#' This method support target set with duplications by transforming it to a 
#' scored ranked target list and conducting enrichment analysis by modified 
#' GSEA method.
#' 
#' The original GSEA method proposed by Subramanian et at., 2005 uses 
#' predefined gene sets \emph{S} given by the GO or KEGG annotations. The goal 
#' is to determine whether the genes in \emph{S} are randomly distributed 
#' throughout a ranked test gene list \emph{L} (e.g. all genes ranked by 
#' log2 fold changes) or enriched at the top or bottom. This is expressed by an 
#' Enrichment Score (\emph{ES}) reflecting the degree to which a set \emph{S} 
#' is overrepresented at the extremes of \emph{L}. 
#' 
#' For TSEA, the query is a target set where duplicated entries are maintained. 
#' To perform GSEA, this target set with duplications was transformed to a 
#' scored ranked target list \emph{L_tar} of all targets provided by the 
#' corresponding annotation system. For each target in the query target set, 
#' its frequency is divided by the number of targets in the target set, 
#' which is the weight of that target. 
#' For targets present in the annotation system but absent in the 
#' target set, their scores are set to 0. Thus, every target in the annotation 
#' system will be assigned a score and then sorted decreasingly to obtain
#' \emph{L_tar}.
#' 
#' In case of TSEA, the original GSEA method cannot be used directly since a 
#' large portion of zeros exists in \emph{L_tar}. If the scores of the genes in 
#' set \emph{S} are all zeros, \emph{N_R} (sum of scores of genes in set 
#' \emph{S}) will be zero, which cannot be used as the denominator. 
#' In this case, \emph{ES} is set to -1. If only some genes in set \emph{S}
#' have scores of zeros then \emph{N_R} is set to a larger number to decrease 
#' the weight of the genes in \emph{S} that have non-zero scores.
#'  
#' The reason for this modification is that if only one gene in gene set 
#' \emph{S} has a non-zero score and this gene ranks high in \emph{L_tar}, 
#' the weight of this gene will be 1 resulting in an \emph{ES(S)} close to 1. 
#' Thus, the original GSEA method will score the gene set \emph{S} as 
#' significantly enriched. However, this is undesirable because in this 
#' example only one gene is shared among the target set and the gene set 
#' \emph{S}. Therefore, giving small weights to genes in \emph{S} that
#' have zero scores could decrease the weight of the genes in \emph{S} that 
#' have non-zero scores, thereby decrease the false positive rate. 
#' @title mGSEA method for TSEA
#' @param drugs query drug set used to do TSEA.
#' Can be top ranking drugs in GESS result. 
#' @param type onoe of `GO` or `KEGG`
#' @param ont if type is `GO`, set ontology, one of `BP`,`MF`,`CC` or `ALL`
#' @param nPerm permutation numbers used to calculate p value
#' @param exponent weight of each step
#' @param pAdjustMethod p value adjustment method,
#' one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param pvalueCutoff p value cutoff
#' @param minGSSize minimum size of each gene set in annotation system
#' @param maxGSSize maximum size of each gene set in annotation system
#' @param verbose TRUE or FALSE, print message or not
#' @return \code{\link{feaResult}} object, 
#' represents enriched functional categories.
#' @seealso \code{\link{feaResult}}, \code{\link{fea}}
#' @references Subramanian et at., 2005,
#' \url{https://www.pnas.org/content/102/43/15545}
#' @examples 
#' data(drugs)
#' ## GO annotation system
#' mgsea_res <- tsea_mGSEA(drugs=drugs, type="GO", ont="MF", exponent=1, 
#'                         nPerm=1000, pvalueCutoff=1, minGSSize=5)
#' result(mgsea_res)
#  ## KEGG annotation system
#' mgsea_k_res <- tsea_mGSEA(drugs=drugs, type="KEGG", exponent=1, 
#'                           nPerm=1000, pvalueCutoff=1, minGSSize=2)
#' result(mgsea_k_res)
#' @export
tsea_mGSEA <- function(drugs, 
                      type="GO", ont="MF", 
                      nPerm=1000, exponent=1, 
                      pAdjustMethod = "BH", pvalueCutoff = 0.05,
                      minGSSize = 2, maxGSSize = 500, verbose=FALSE){
  drugs <- unique(tolower(drugs))
  targets <- get_targets(drugs, database = "all")
  gnset <- na.omit(unlist(lapply(targets$t_gn_sym, function(i) 
    unlist(strsplit(as.character(i), split = "; ")))))
  # give scores to gnset
  tar_tab <- table(gnset)
  tar_dup <- as.numeric(tar_tab); names(tar_dup) <- names(tar_tab)
  tar_weight <- sort(tar_dup/sum(tar_dup), decreasing = TRUE)
  
  if(type=="GO"){
    # Get universe genes in GO annotation system
    ext_path <- system.file("extdata", package="signatureSearch")
    univ_tar_path <- paste0(ext_path,"/univ_genes_in_GO.txt")
    universe <- readLines(univ_tar_path)
    tar_diff <- setdiff(universe, gnset)
    tar_diff_weight <- rep(0, length(tar_diff))
    names(tar_diff_weight)=tar_diff
    tar_total_weight = c(tar_weight, tar_diff_weight)
    gsego <- gseGO2(geneList = tar_total_weight, OrgDb = org.Hs.eg.db, 
                    ont = ont, keyType = "SYMBOL", nPerm = nPerm, 
                    minGSSize = minGSSize, maxGSSize = maxGSSize, 
                    exponent = exponent, nproc=1, verbose=verbose, 
                    pvalueCutoff = pvalueCutoff, pAdjustMethod=pAdjustMethod)
    if(is.null(gsego))
        return(NULL)
    gsego@drugs = drugs
    return(gsego)
  }
  
  if(type=="KEGG"){
    # Get universe genes in KEGG annotation system
    ext_path <- system.file("extdata", package="signatureSearch")
    univ_tar_path <- paste0(ext_path,"/univ_genes_in_KEGG.txt")
    universe <- readLines(univ_tar_path)
    
    # convert gnset symbol to entrez
    gnset_entrez <- suppressMessages(
      AnnotationDbi::select(org.Hs.eg.db, keys = gnset, 
                            keytype = "SYMBOL", columns = "ENTREZID"))
    gnset_entrez2 <- as.character(na.omit(gnset_entrez$ENTREZID))
    # give scores to gnset_entrez2
    tar_tab <- table(gnset_entrez2)
    tar_dup <- as.numeric(tar_tab); names(tar_dup) <- names(tar_tab)
    tar_weight <- sort(tar_dup/sum(tar_dup), decreasing = TRUE)
    
    tar_diff <- setdiff(universe, gnset_entrez2)
    tar_diff_weight <- rep(0, length(tar_diff))
    names(tar_diff_weight)=tar_diff
    tar_total_weight = c(tar_weight, tar_diff_weight)
    
    gsekk <- gseKEGG2(geneList=tar_total_weight, organism='hsa', 
                      keyType='kegg', nPerm = nPerm, verbose = verbose, 
                      minGSSize = minGSSize, maxGSSize=maxGSSize, 
                      pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod)
    if(is.null(gsekk))
        return(NULL)
    gsekk@drugs = drugs
    return(gsekk)
  }
}
