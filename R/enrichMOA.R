##' Given a vector of gene identifiers, this function returns MOA category 
##' enrichment results based on a hypergeometric test with duplication support 
##' in the test set. The universe for the test is set to the unique genes 
##' encoding the target proteins present in the MOA annotation system from the 
##' ChEMBL database. 
##' 
##' @title MOA Category Enrichment with Hypergeometric Test
##' @param gene a vector of entrez gene ids (here the test set)
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", 
##' "bonferroni", "BH", "BY", "fdr", "none"
##' @param qvalueCutoff qvalue cutoff
##' @return A \code{feaResult} instance.
##' @seealso \code{\link{feaResult-class}}
##' @examples 
##' data(geneList, package="DOSE")
##' emoa <- enrichMOA(gene = names(geneList)[seq(3)])
##' head(emoa)
##' @export
enrichMOA <- function(gene,
                      pvalueCutoff=0.05,
                      pAdjustMethod="BH",
                      qvalueCutoff = 0.2) {
    data("chembl_moa_list", envir = environment())
    MOA_DATA_chembl <- get_MOA_data(chembl_moa_list, keytype="entrez")
    # get all the gene entrez ids in the MOA annotation system as universe
    ext2path <- get("EXTID2PATHID", envir = MOA_DATA_chembl)
    universe <- names(ext2path)
    
    res <- enricher_internal(gene,
                           pvalueCutoff=pvalueCutoff,
                           pAdjustMethod=pAdjustMethod,
                           universe = universe,
                           qvalueCutoff = qvalueCutoff,
                           minGSSize = 0,
                           maxGSSize = 100,
                           USER_DATA = MOA_DATA_chembl)
    
    if (is.null(res))
        return(res)
    og(res) <- "Homo sapiens"
    ont(res) <- "MOA"
    return(res)
}
