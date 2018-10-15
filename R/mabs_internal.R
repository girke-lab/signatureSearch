
##' generic function for meanAbs analysis
##'
##' @title mabs_internal
##' @param geneList order ranked geneList
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of each geneSet for analyzing
##' @param pvalueCutoff p value Cutoff
##' @param pAdjustMethod p value adjustment method
##' @param USER_DATA annotation data
##' @importFrom stats p.adjust
##' @importFrom qvalue qvalue
##' @return mabsResult object
##' @author Yuzhu Duan
mabs_internal <- function(geneList,
                          nPerm,
                          minGSSize,
                          maxGSSize,
                          pvalueCutoff,
                          pAdjustMethod,
                          USER_DATA) {
    geneList <- sort(geneList, decreasing = TRUE)
    geneSets <- getGeneSet(USER_DATA)
    selected.gs <- geneSet_filter(geneSets, geneList, minGSSize, maxGSSize)

    if (is.null(selected.gs))
        return(NULL)
  
    message("calculating observed meanAbs scores...")

    observedScore <- sapply(selected.gs, function(gs)
        mabs_score(geneSet=gs, geneList=geneList)
        )
    
    # exclude gene sets that have no intersect with drug targets
    selected.gs <- selected.gs[observedScore != 0]
    observedScore <- observedScore[observedScore != 0]
    
    message("calculating permutation scores...")
    
    permScores <- lapply(1:nPerm, function(i) {
        perm.mabs_score(geneList, selected.gs)
    })

    permScores <- do.call("cbind", permScores)

    rownames(permScores) <- names(selected.gs)

    median <- apply(permScores, 1, median)
    sd <- apply(permScores, 1, sd)

    Nmabs <- (observedScore - median)/sd
    
    message("calculating p values...")
    
    pvals <- sapply(seq_along(observedScore), function(i) {
        sum(permScores[i,]>observedScore[i]) / nPerm
    })
    
    p.adj <- p.adjust(pvals, method=pAdjustMethod)
    qvalues <- calculate_qvalue(pvals)

    gs.name <- names(selected.gs)
    Description <- TERM2NAME(gs.name, USER_DATA)

    geneID <- sapply(selected.gs, function(x) paste0(intersect(x, names(geneList)[geneList>0]), collapse = "/"))
      
    params <- list(pvalueCutoff = pvalueCutoff,
                   nPerm = nPerm,
                   pAdjustMethod = pAdjustMethod,
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize
                   )


    res <- data.frame(
        ID = as.character(gs.name),
        Description = Description,
        setSize = sapply(selected.gs, length),
        mabs = observedScore,
        Nmabs = Nmabs,
        pvalue = pvals,
        p.adjust = p.adj,
        qvalues = qvalues,
        geneID = geneID,
        stringsAsFactors = FALSE
    )

    res <- res[!is.na(res$pvalue),]
    res <- res[ res$pvalue <= pvalueCutoff, ]
    res <- res[ res$p.adjust <= pvalueCutoff, ]
    idx <- order(res$Nmabs, decreasing = TRUE)
    res <- res[idx, ]
    row.names(res) <- res$ID
    if (nrow(res) == 0) {
        message("no term enriched under specific pvalueCutoff...")
        res = new("mabsResult",
                result     = res,
                geneSets   = geneSets,
                geneList   = geneList,
                params     = params,
                readable   = TRUE
                )
    }
    message("done...")

    res = new("feaResult",
        result     = as_tibble(res),
        refSets   = geneSets,
        targets   = geneList,
        universe = names(geneList)
        )
    res@organism <- "UNKNOWN"
    return(res)
}

mabs_score <- function(geneList, geneSet) {
    mabs <- mean(geneList[geneSet])
    return(mabs)
}

perm.geneList <- function(geneList) {
    ## perm.idx <- sample(seq_along(geneList), length(geneList), replace=FALSE)
    perm.idx <- sample.int(length(geneList))
    perm.geneList <- geneList
    names(perm.geneList) <- names(geneList)[perm.idx]
    return(perm.geneList)
}

perm.mabs_score <- function(geneList, geneSets) {
    geneList <- perm.geneList(geneList)
    res <- sapply(1:length(geneSets), function(i)
                  mabs_score(geneSet=geneSets[[i]],
                             geneList=geneList)
                  )
    return(res)
}


geneSet_filter <- function(geneSets, geneList, minGSSize, maxGSSize) {
    geneSets <- sapply(geneSets, intersect, names(geneList))

    gs.idx <- get_geneSet_index(geneSets, minGSSize, maxGSSize)
    nGeneSet <- sum(gs.idx)

    if ( nGeneSet == 0 ) {
        msg <- paste0("No gene set have size between [", minGSSize, ", ", maxGSSize, "]...")
        message(msg)
        message("--> return NULL...")
        return(NULL)
    }
    geneSets[gs.idx]
}

