##' Given a vector of genes, this function will return the enriched GO
##' categories after FDR control.
##' The hypergeometric test is adjusted to support gene set with duplications
##' 
##' @title GO enrichment analysis via hypergeometric test
##' @param gene a vector of entrez gene id or gene SYMBOL.
##' @param OrgDb OrgDb
##' @param keytype keytype of input gene
##' @param ont One of "MF", "BP", and "CC" subontologies.
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
  
  if (ont == "ALL" && !pool) {
    lres <- lapply(c("BP", "CC", "MF"), function(ont)
      suppressMessages(enrichGO2(gene, OrgDb, keytype, ont,
                                pvalueCutoff, pAdjustMethod, universe,
                                qvalueCutoff, minGSSize, maxGSSize
      ))
    )
    
    lres <- lres[!vapply(lres, is.null, logical(1))]
    if (length(lres) == 0)
      return(NULL)
    
    df <- do.call('rbind', lapply(lres, as.data.frame))
    refSets <- lres[[1]]@refSets
    if (length(lres) > 1) {
      for (i in 2:length(lres)) {
        refSets <- append(refSets, lres[[i]]@refSets)
      }
    }
    res <- lres[[1]]
    res@result <- df
    res@refSets <- refSets
  } else {
    res <- enricher_internal(gene,
                             pvalueCutoff=pvalueCutoff,
                             pAdjustMethod=pAdjustMethod,
                             universe = universe,
                             qvalueCutoff = qvalueCutoff,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             USER_DATA = GO_DATA
    )
    
    if (is.null(res))
      return(res)
  }
  res@organism <- DOSE:::get_organism(OrgDb)
  res@ontology <- ont
  
  if (ont == "ALL") {
    res <- clusterProfiler:::add_GO_Ontology(res, GO_DATA)
  }
  return(res)
}
