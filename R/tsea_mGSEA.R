#' The Modified Gene Set Enrichment Analysis (mGSEA) method supports
#' target set with duplications by transforming it to a score ranked target 
#' list and conducting the enrichment analysis with GSEA algorithm with
#' modification.
#' 
#' The original GSEA method proposed by Subramanian et at., 2005 uses 
#' predefined gene sets \emph{S} defined by functional annotation systems 
#' such as GO and KEGG. The goal is to determine whether the genes in \emph{S}
#' are randomly distributed throughout a ranked test gene list \emph{L} 
#' (e.g. all genes ranked by log2 fold changes) or enriched at the top or 
#' bottom of the test lit. This is expressed by an 
#' Enrichment Score (\emph{ES}) reflecting the degree to which a set \emph{S} 
#' is overrepresented at the extremes of \emph{L}. 
#' 
#' For TSEA, the query is a target set where duplicated entries need to be
#' maintained. To perform GSEA with duplication support, here referred to as 
#' mGSEA, the target set is transformed to a score ranked target list 
#' \emph{L_tar} of all targets provided by the 
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
#' To favor truly enriched GO terms and KEGG pathways (gene set \emph{S}) at 
#' the top of \emph{L_tar}, only gene sets with positive \emph{ES} are selected.
#' @section Column description:
#' Description of the columns in the result table specific to the GSEA 
#' algorithm:
#' \itemize{
#'     \item enrichmentScore: ES from the GSEA algorithm 
#'     (Subramanian et al., 2005). The score is calculated by walking down the 
#'     gene list L, increasing a running-sum statistic when we encounter a gene 
#'     in S and decreasing when it is not. The magnitude of the increment 
#'     depends on the gene scores. The ES is the maximum deviation from zero 
#'     encountered in the random walk. It corresponds to a weighted 
#'     Kolmogorov-Smirnov-like statistic.
#'     \item NES: Normalized enrichment score. The positive and negative 
#'     enrichment scores are normalized separately by permutating the 
#'     composition of the gene list L nPerm times, and dividing the enrichment 
#'     score by the mean of the permutation ES with the same sign.
#'     \item pvalue:  The nominal p-value of the ES is calculated using a 
#'     permutation test. Specifically, the composition of the gene list L is 
#'     permuted and the ES of the gene set is recomputed for the permutated 
#'     data generating a null distribution for the ES. The p-value of the 
#'     observed ES is then calculated relative to this null distribution.
#'     \item leadingEdge: Genes in the gene set S (functional category) that 
#'     appear in the ranked list L at, or before, the point where the running 
#'     sum reaches its maximum deviation from zero. It can be interpreted as 
#'     the core of a gene set that accounts for the enrichment signal.
#'     \item ledge_rank: Ranks of genes in 'leadingEdge' in gene list L.
#' }
#' Description of the other columns are available at the 'result' slot of the
#' \code{\link{feaResult}} object.
#' 
#' @title mGSEA Enrichment Method
#' @param drugs character vector, query drug set used for functional enrichment.
#' Can be top ranking drugs in the GESS result. 
#' @param type one of `GO` or `KEGG`
#' @param ont character(1). If type is `GO`, set ontology as one of `BP`,`MF`,
#' `CC` or `ALL`. If type is 'KEGG', it is ignored.
#' @param nPerm integer, permutation numbers used to calculate p-value
#' @param exponent integer, exponent in GSEA algorithm defines the weight of 
#' genes in gene set \emph{S}.
#' @param pAdjustMethod p-value adjustment method, 
#' one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'
#' @param pvalueCutoff double, p-value cutoff
#' @param minGSSize integer, minimum size of each gene set in annotation system
#' @param maxGSSize integer, maximum size of each gene set in annotation system
#' @param verbose TRUE or FALSE, print message or not
#' @param dt_anno drug-target annotation resource. one of 'DrugBank', 'CLUE', 
#' 'STITCH' or 'all'. If 'dt_anno' is 'all', the targets from DrugBank, CLUE 
#' and STITCH databases will be combined. It is recommended to set the 'dt_anno'
#' as 'all' since it will get the most complete target set of as many drugs
#' as possible. Users could also choose individual annotation resource if 
#' wanted, but should be aware that if the chosen drug-target annotation
#' resource contains limited drugs (such as CLUE), many query drugs will not
#' get targets and the target set may not be complete, which could affect
#' the enrichment result.  
#' @return \code{\link{feaResult}} object, the result table contains the
#' enriched functional categories (e.g. GO terms or KEGG pathways) ranked by 
#' the corresponding enrichment statistic.
#' @seealso \code{\link{feaResult}}, \code{\link{fea}}
#' @references 
#' GSEA algorithm: 
#' Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., 
#' Gillette, M. A., … Mesirov, J. P. (2005). Gene set enrichment analysis: a 
#' knowledge-based approach for interpreting genome-wide expression profiles. 
#' Proceedings of the National Academy of Sciences of the United States of
#' America, 102(43), 15545–15550. \url{https://doi.org/10.1073/pnas.0506580102}
#' @examples 
#' data(drugs)
#' ## GO annotation system
#' #mgsea_res <- tsea_mGSEA(drugs=drugs, type="GO", ont="MF", exponent=1, 
#' #                        nPerm=1000, pvalueCutoff=1, minGSSize=5)
#' #result(mgsea_res)
#  ## KEGG annotation system
#' mgsea_k_res <- tsea_mGSEA(drugs=drugs, type="KEGG", exponent=1, 
#'                           nPerm=1000, pvalueCutoff=1, minGSSize=2)
#' result(mgsea_k_res)
#' @export
tsea_mGSEA <- function(drugs, 
                      type="GO", ont="MF", 
                      nPerm=1000, exponent=1, 
                      pAdjustMethod="BH", pvalueCutoff=0.05,
                      minGSSize=2, maxGSSize=500, verbose=FALSE,
                      dt_anno="all"){
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
