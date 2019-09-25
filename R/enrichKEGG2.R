##' Given a vector of gene identifiers, this function returns KEGG pathway 
##' enrichment results based on a hypergeometric test with duplication support 
##' in the test set.
##' 
##' @title KEGG Pathway Enrichment with Hypergeometric Test
##' @param gene a vector of entrez gene ids (here the test set)
##' @param organism supported organism are listed in 
##' http://www.genome.jp/kegg/catalog/org_list.html
##' @param keyType one of "kegg", 'ncbi-geneid', 'ncbi-proteinid' or 'uniprot'
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", 
##' "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by 
##' ontology term for testing.
##' @param maxGSSize maximal size of genes annotated for testing
##' @param qvalueCutoff qvalue cutoff
##' @return A \code{feaResult} instance.
##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importMethodsFrom AnnotationDbi mget
##' @importClassesFrom methods data.frame
##' @examples 
##' # Method supports duplicated elements in "gene", which should be entrez ids
##' gene <- c(rep("4312",4), rep("8318",2), "991", "10874")
##' #data(geneList, package="DOSE")
##' #kk <- enrichKEGG2(gene = gene, universe=names(geneList))
##' #head(kk)
##' @export
enrichKEGG2 <- function(gene,
                       organism          = "hsa",
                       keyType           = "kegg",
                       pvalueCutoff      = 0.05,
                       pAdjustMethod     = "BH",
                       universe,
                       minGSSize         = 5,
                       maxGSSize         = 500,
                       qvalueCutoff      = 0.2) {
  
  species <- organismMapper(organism)
  KEGG_DATA <- prepare_KEGG(species, "KEGG", keyType)
  res <- enricher_internal(gene,
                           pvalueCutoff  = pvalueCutoff,
                           pAdjustMethod = pAdjustMethod,
                           universe      = universe,
                           minGSSize     = minGSSize,
                           maxGSSize     = maxGSSize,
                           qvalueCutoff  = qvalueCutoff,
                           USER_DATA = KEGG_DATA)
  if (is.null(res))
    return(res)
  ont(res) <- "KEGG"
  og(res) <- species
  return(res)
}
