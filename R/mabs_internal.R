## generic function for meanAbs analysis
mabs_internal <- function(geneList,
                          nPerm,
                          minGSSize,
                          maxGSSize,
                          pvalueCutoff,
                          pAdjustMethod,
                          USER_DATA) {
    geneList <- sort(geneList, decreasing = TRUE)
    geneSets <- get("PATHID2EXTID", envir = USER_DATA)
    selected.gs <- geneSet_filter(geneSets, geneList, minGSSize, maxGSSize)

    if (is.null(selected.gs))
        return(NULL)
  
    #message("calculating observed meanAbs scores...")

    observedScore <- vapply(selected.gs, function(gs)
        mabs_score(geneSet=gs, geneList=geneList),
        FUN.VALUE = numeric(1))
    
    # exclude gene sets that have no intersect with drug targets
    selected.gs <- selected.gs[observedScore != 0]
    observedScore <- observedScore[observedScore != 0]
    
    #message("calculating permutation scores...")
    
    # get 1000 permutation matrix of geneList
    perm_mat <- replicate(nPerm, sample(geneList))
    permScores <- vapply(selected.gs, function(gs){
      colMeans(abs(perm_mat[gs,]))},
      FUN.VALUE = numeric(ncol(perm_mat)))

    median <- apply(permScores, 2, median)
    sd <- apply(permScores, 2, sd)

    Nmabs <- (observedScore - median)/sd
    
    #message("calculating p values...")
    
    pvals <- vapply(seq_along(observedScore), function(i) {
        sum(permScores[,i]>observedScore[i]) / nPerm
    }, FUN.VALUE = numeric(1))
    
    p.adj <- p.adjust(pvals, method=pAdjustMethod)
    qvalues <- calculate_qvalue(pvals)

    gs.name <- names(selected.gs)
    Description <- TERM2NAME(gs.name, USER_DATA)

    geneID <- vapply(selected.gs, function(x) 
      paste0(intersect(x, names(geneList)[geneList>0]), collapse = "/"),
      FUN.VALUE = character(1))

    res <- data.frame(
        ID = as.character(gs.name),
        Description = Description,
        setSize = vapply(selected.gs, length, FUN.VALUE = integer(1)),
        mabs = observedScore,
        Nmabs = Nmabs,
        pvalue = pvals,
        p.adjust = p.adj,
        qvalues = qvalues,
        itemID = geneID,
        stringsAsFactors = FALSE
    )

    res <- res[!is.na(res$pvalue),]
    res <- res[ res$pvalue <= pvalueCutoff, ]
    res <- res[ res$p.adjust <= pvalueCutoff, ]
    # order by mabs
    idx <- order(-res$mabs)
    res <- res[idx, ]
    row.names(res) <- res$ID
    if (nrow(res) == 0) {
        message("no term enriched under specific pvalueCutoff...")
        return(NULL)
    }
    #message("done...")

    res = new("feaResult",
        result     = as_tibble(res),
        #refSets   = geneSets,
        targets   = geneList
        #universe = names(geneList)
        )
    res@organism <- "UNKNOWN"
    return(res)
}

mabs_score <- function(geneList, geneSet) {
    mabs <- mean(abs(geneList[geneSet]))
    return(mabs)
}

