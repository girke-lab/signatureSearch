
GSEA_fgsea2 <- function(geneList,
                       exponent,
                       nPerm,
                       minGSSize,
                       maxGSSize,
                       pvalueCutoff,
                       pAdjustMethod,
                       verbose,
                       nproc=1,
                       seed=FALSE,
                       USER_DATA) {

    if(verbose)
        message("preparing geneSet collections...")

    geneSets <- getGeneSet(USER_DATA)
    
    if(verbose)
      message("excluding gene sets that have no intersect with drug targets")
    
    logic <- sapply(geneSets, function(x) ifelse(length(intersect(x, names(geneList)[geneList!=0]))==0, FALSE, TRUE))
    geneSets = geneSets[logic]
    
    if(verbose)
      message("Filtering gene sets for genes not in universe and excluding gene sets that are less than 'minGSSize' or greater than 'maxGSSize'")
    geneSets <- sapply(geneSets, function(x) intersect(x, names(geneList)))
    logic <- sapply(geneSets, function(x) length(x)<=maxGSSize & length(x)>=minGSSize)
    geneSets = geneSets[logic]

    if(verbose)
        message("GSEA analysis...")

    tmp_res <- fgsea2(pathways=geneSets,
                 stats=geneList,
                 nperm=nPerm,
                 minSize=minGSSize,
                 maxSize=maxGSSize,
                 gseaParam=exponent,
                 nproc = nproc)

    p.adj <- p.adjust(tmp_res$pval, method=pAdjustMethod)
    qvalues <- calculate_qvalue(tmp_res$pval)

    Description <- TERM2NAME(tmp_res$pathway, USER_DATA)

    ledge <- sapply(tmp_res$leadingEdge, paste0, collapse='/')
    ledge_rank <- sapply(tmp_res$ledge_rank, paste0, collapse='/')
    message("ledge_rank included")
    res <- data.frame(
        ID = as.character(tmp_res$pathway),
        Description = Description,
        setSize = tmp_res$size,
        enrichmentScore = tmp_res$ES,
        NES = tmp_res$NES,
        pvalue = tmp_res$pval,
        p.adjust = p.adj,
        qvalues = qvalues,
        leadingEdge = ledge,
        ledge_rank = ledge_rank,
        stringsAsFactors = FALSE
    )

    res <- res[!is.na(res$pvalue),]
    res <- res[ res$pvalue <= pvalueCutoff, ]
    res <- res[ res$p.adjust <= pvalueCutoff, ]
    idx <- order(res$enrichmentScore, decreasing = TRUE)
    res <- res[idx, ]

    if (nrow(res) == 0) {
        message("no term enriched under specific pvalueCutoff...")
        return(
            new("feaResult",
                result     = res,
                refSets   = geneSets,
                targets   = geneList,
                universe  = names(geneList)
                )
        )
    }

    row.names(res) <- res$ID
    
    if (verbose)
        message("done...")

    new("feaResult",
        result    = as_tibble(res),
        refSets   = geneSets,
        targets   = geneList,
        universe  = names(geneList)
        )
}

##' generic function for gene set enrichment analysis (drug targets score list as geneList,
##' GSEA method is modified to accept geneList with large portion of zeros)
##'
##' @title GSEA_internal2
##' @param geneList order ranked geneList
##' @param exponent weight of each step
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of each geneSet for analyzing
##' @param pvalueCutoff p value Cutoff
##' @param pAdjustMethod p value adjustment method
##' @param verbose print message or not
##' @param seed set seed inside the function to make result reproducible. FALSE by default.
##' @param USER_DATA annotation data
##' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 1)
##' @param by one of 'fgsea' or 'DOSE'
##' @return feaResult object
GSEA_internal2 <- function(geneList,
                 exponent,
                 nPerm,
                 minGSSize,
                 maxGSSize,
                 pvalueCutoff,
                 pAdjustMethod,
                 verbose,
                 seed=FALSE,
                 USER_DATA,
                 nproc=1,
                 by="fgsea") {

    by <- match.arg(by, c("fgsea", "DOSE"))
    if (is.unsorted(-geneList))
        stop("geneList should be a decreasing sorted vector...")
    if (by == 'fgsea') {
        .GSEA <- GSEA_fgsea2
    } 
    res <- .GSEA(geneList          = geneList,
                 exponent          = exponent,
                 nproc             = nproc,
                 nPerm             = nPerm,
                 minGSSize         = minGSSize,
                 maxGSSize         = maxGSSize,
                 pvalueCutoff      = pvalueCutoff,
                 pAdjustMethod     = pAdjustMethod,
                 verbose           = verbose,
                 seed              = seed,
                 USER_DATA         = USER_DATA)
    res@organism <- "UNKNOWN"
    return(res)
}