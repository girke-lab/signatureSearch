#' meanAbs method for TSEA
#' 
#' This method support target set with duplications by transforming it to a 
#' scored ranked target list and calculating mean absolute scores of genes
#' in gene set \emph{S}
#' 
#' The input for the \emph{mabs} method is \emph{L_tar}, the same as for 
#' \code{m_GSEA}. The \emph{meanAbs} statistic (\emph{mabs(S)}) of a gene set 
#' \emph{S} is calculated as mean absolute scores of the genes in \emph{S}. 
#' In order to adjust for size variations in gene set \emph{S}, 1000 random 
#' permutations of \emph{L_tar} are performed to determine \emph{mabs(S,pi)}. 
#' Subsequently, \emph{mabs(S)} is normalized by subtracting the median of the 
#' \emph{mabs(S,pi)} and then dividing by the standard deviation of 
#' \emph{mabs(S,pi)} yielding the normalized scores \emph{Nmabs(S)}. Finally, 
#' the portion of \emph{mabs(S,pi)} that is greater than \emph{mabs(S)} is 
#' used as nominal p value (Fang et al., 2012). 
#' The resulting nominal p values are adjusted for 
#' multiple hypothesis testing using the Benjamini-Hochberg method.
#' 
#' Description of the columns in the result table specific to the MeanAbs 
#' algorithm:
#' \itemize{
#'     \item mabs: Given a scored ranked gene list L, mabs(S) represents
#'     the mean absolute scores of the genes in set S. 
#'     \item Nmabs: mabs(S) normalized
#' }
#' @param drugs query drug set used to do TSEA. 
#' Can be top ranking drugs in GESS result. 
#' @param type can be `GO` or `KEGG`
#' @param ont if type is `GO`, set ontology, can be `BP`,`MF`,`CC` or `ALL`
#' @param nPerm permutation numbers used to calculate p value
#' @param pAdjustMethod p value adjustment method,
#' one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param pvalueCutoff p value cutoff
#' @param minGSSize minimum size of each gene set in annotation system
#' @param maxGSSize maximum size of each gene set in annotation system
#' @return \code{\link{feaResult}} object, 
#' represents enriched functional categories.
#' @seealso \code{\link{feaResult}}, \code{\link{fea}}, \code{\link{tsea_mGSEA}}
#' @references Fang, Z., Tian, W., & Ji, H. (2012). A network-based 
#' gene-weighting approach for pathway analysis. Cell Research, 22(3), 
#' 565–580. \url{https://doi.org/10.1038/cr.2011.149}
#' @examples 
#' data(drugs)
#' ## GO annotation system
#' #mabs_res <- tsea_mabs(drugs=drugs, type="GO", ont="MF", nPerm=1000, 
#' #                      pvalueCutoff=0.05, minGSSize=5)
#' #result(mabs_res)
#' ## KEGG annotation system
#' #mabs_k_res <- tsea_mabs(drugs=drugs, type="KEGG", nPerm=1000, 
#' #                        pvalueCutoff=0.05, minGSSize=5)
#' #result(mabs_k_res)
#' @export
#' 
tsea_mabs <- function(drugs, 
                      type="GO", ont="MF", 
                      nPerm=1000,  
                      pAdjustMethod = "BH", pvalueCutoff = 0.05,
                      minGSSize = 5, maxGSSize = 500){
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
    mabsgo <- mabsGO(geneList = tar_total_weight, OrgDb = org.Hs.eg.db, 
                     ont = ont, keyType = "SYMBOL", nPerm = nPerm, 
                     minGSSize = minGSSize, maxGSSize = maxGSSize, 
                     pvalueCutoff = pvalueCutoff, pAdjustMethod=pAdjustMethod)
    mabsgo@drugs = drugs
    return(mabsgo)
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
    tar_diff_weight <- rep(0, length(tar_diff)); names(tar_diff_weight)=tar_diff
    tar_total_weight = c(tar_weight, tar_diff_weight)

    mabskk <- mabsKEGG(geneList=tar_total_weight, organism='hsa', 
                       keyType='kegg', nPerm = nPerm, 
                       minGSSize = minGSSize, maxGSSize=maxGSSize, 
                       pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod)
    if(is.null(mabskk))
      return(NULL)
    mabskk@drugs = drugs
    return(mabskk)
  }
}