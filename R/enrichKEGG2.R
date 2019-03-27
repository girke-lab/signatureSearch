##' Given a vector of genes, this function will return the enriched KEGG
##' pathways after FDR control.
##' The hypergeometric test is adjusted to support gene set with duplications
##' 
##' @title KEGG pathways enrichment analysis via hypergeometric test
##' @param gene a vector of entrez gene id.
##' @param organism supported organism listed in 
##' \url{http://www.genome.jp/kegg/catalog/org_list.html}
##' @param keyType one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
##' @param pvalueCutoff Cutoff value of pvalue.
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
  
  species <- clusterProfiler:::organismMapper(organism)
  KEGG_DATA <- clusterProfiler:::prepare_KEGG(species, "KEGG", keyType)
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
  res@ontology <- "KEGG"
  res@organism <- species
  return(res)
}
