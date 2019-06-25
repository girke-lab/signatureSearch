##' Given a vector of gene entrez ids, this function will return the enriched 
##' MOA categories after FDR control. The universe for the hypergeometric test
##' are all the genes in the MOA annotation obtained from ChEMBL database.
##' 
##' @title MOA enrichment analysis via hypergeometric test
##' @param gene a vector of gene entrez ids.
##' @param pvalueCutoff Cutoff value of pvalue.
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
    chembl_moa_list <- readRDS(system.file("extdata", 
                                           "ChEMBL_moa2entrezid_list.rds", 
                                    package="signatureSearch"))
    MOA_DATA_chembl <- get_MOA_data(chembl_moa_list, keytype="entrez")
    # get all the gene entrez ids in the MOA annotation system as universe
    ext2path <- get("EXTID2PATHID", envir = MOA_DATA_chembl)
    universe = names(ext2path)
    
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
    res@organism <- "Homo sapiens"
    res@ontology <- "MOA"
    return(res)
}
