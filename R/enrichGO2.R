##' Given a vector of genes, this function will return the enriched GO
##' categories after FDR control.
##' The hypergeometric test is adjusted to support gene set with duplications
##' 
##' @title GO enrichment analysis via hypergeometric test
##' @param gene a vector of entrez gene id or gene SYMBOL.
##' @param OrgDb OrgDb
##' @param keytype keytype of input gene
##' @param ont One of "MF", "BP", "CC" or "ALL"
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", 
##' "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param qvalueCutoff qvalue cutoff
##' @param minGSSize minimal size of genes annotated by Ontology term 
##' @param maxGSSize maximal size of genes annotated for testing
##' @param pool If ont='ALL', whether pool 3 GO sub-ontologies
##' @return A \code{feaResult} instance.
##' @seealso \code{\link{feaResult-class}}
##' @examples 
##' # The method supports duplicated elements in 'gene', 
##' # which should be SYMBOL id for GO enrichment.
##' gene = c(rep("HDAC1",4), rep("HDAC3",2), "SOX8", "KLK14")
##' library(org.Hs.eg.db)
##' data(targetList)
##' ego <- enrichGO2(gene = gene, OrgDb=org.Hs.eg.db, ont="MF",
##'                  universe=names(targetList))
##' head(ego)
##' @export
enrichGO2 <- function(gene,
                     OrgDb,
                     keytype = "SYMBOL",
                     ont="MF",
                     pvalueCutoff=0.05,
                     pAdjustMethod="BH",
                     universe,
                     qvalueCutoff = 0.2,
                     minGSSize = 5,
                     maxGSSize = 500,
                     pool=FALSE) {
  
  ont %<>% toupper
  ont <- match.arg(ont, c("BP", "CC", "MF", "ALL"))
  # GO_DATA <- clusterProfiler:::get_GO_data(OrgDb, ont, keytype)
  # download GO_DATA.rds and save it to cache to save time
  fl <- download_data_file(url=
        "http://biocluster.ucr.edu/~yduan004/signatureSearch_data/GO_DATA.rds",
        rname="GO_DATA")
  GO_DATA <- readRDS(fl)
  if (missing(universe))
    universe <- NULL
  
  res <- enricher_internal(gene,
                           pvalueCutoff=pvalueCutoff,
                           pAdjustMethod=pAdjustMethod,
                           universe = universe,
                           qvalueCutoff = qvalueCutoff,
                           minGSSize = minGSSize,
                           maxGSSize = maxGSSize,
                           USER_DATA = GO_DATA)
    
  if (is.null(res))
    return(res)
  # Add and select ontology in res
  res <- add_GO_Ontology(res, GO_DATA)
  tmp_df <- res@result
  colnames(tmp_df)[1] = "ont"
  res@result <- tmp_df
  if(ont != "ALL")
    res@result <- as_tibble(res[res$ont == ont, ])
  res@organism <- get_organism(OrgDb)
  res@ontology <- ont
  return(res)
}
