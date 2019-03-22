enrichGO <- function(gene,
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
  # download GO_DATA and save it locally to save time
  ext_path <- system.file("extdata", package="signatureSearch")
  godata_path <- paste0(ext_path,"/GO_DATA.rds")
  if(file.exists(godata_path)){
    GO_DATA <- readRDS(godata_path)
  } else {
    download.file("http://biocluster.ucr.edu/~yduan004/fea/GO_DATA.rds",
    godata_path, quiet = TRUE)
    GO_DATA <- readRDS(godata_path)
  }
  
  if (missing(universe))
    universe <- NULL
  
  if (ont == "ALL" && !pool) {
    lres <- lapply(c("BP", "CC", "MF"), function(ont)
      suppressMessages(enrichGO(gene, OrgDb, keytype, ont,
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
  res@organism <- clusterProfiler:::get_organism(OrgDb)
  res@ontology <- ont
  
  if (ont == "ALL") {
    res <- clusterProfiler:::add_GO_Ontology(res, GO_DATA)
  }
  return(res)
}
