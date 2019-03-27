#' @importFrom fgsea fgsea
GSEA_fgsea <- function(geneList,
                       exponent,
                       nPerm,
                       minGSSize,
                       maxGSSize,
                       pvalueCutoff,
                       pAdjustMethod,
                       verbose,
                       seed=FALSE,
                       USER_DATA) {

    if(verbose)
        message("preparing geneSet collections...")

    geneSets <- get("PATHID2EXTID", envir = USER_DATA)
    DOSE:::check_gene_id(geneList, geneSets)

    if(verbose)
        message("GSEA analysis...")

    tmp_res <- fgsea(pathways=geneSets,
                 stats=geneList,
                 nperm=nPerm,
                 minSize=minGSSize,
                 maxSize=maxGSSize,
                 gseaParam=exponent,
                 nproc = 1)

    p.adj <- p.adjust(tmp_res$pval, method=pAdjustMethod)
    qvalues <- DOSE:::calculate_qvalue(tmp_res$pval)

    Description <- DOSE:::TERM2NAME(tmp_res$pathway, USER_DATA)

    params <- list(pvalueCutoff = pvalueCutoff,
                   nPerm = nPerm,
                   pAdjustMethod = pAdjustMethod,
                   exponent = exponent,
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize
                   )
    ledge <- vapply(tmp_res$leadingEdge, paste0, collapse='/', 
                    FUN.VALUE = character(1))
    ledge_rank <- lapply(tmp_res$leadingEdge, 
                         function(x) match(x, names(geneList)))
    ledge_rank2 <- vapply(ledge_rank, paste, collapse="/", 
                          FUN.VALUE = character(1))
    if(verbose)
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
        ledge_rank = ledge_rank2,
        stringsAsFactors = FALSE
    )

    res <- res[!is.na(res$pvalue),]
    res <- res[ res$pvalue <= pvalueCutoff, ]
    res <- res[ res$p.adjust <= pvalueCutoff, ]
    idx <- order(res$NES, decreasing = TRUE)
    res <- res[idx, ]

    if (nrow(res) == 0) {
        message("No term enriched under specific pvalueCutoff...")
        return(NULL)
    }

    row.names(res) <- res$ID

    if (verbose)
        message("done...")

    new("feaResult",
        result    = as_tibble(res),
        #refSets   = geneSets,
        targets   = geneList
        #universe  = names(geneList)
    )
}

## generic function for gene set enrichment analysis
GSEA_internal <- function(geneList,
                 exponent,
                 nPerm,
                 minGSSize,
                 maxGSSize,
                 pvalueCutoff,
                 pAdjustMethod,
                 verbose=FALSE,
                 seed=FALSE,
                 USER_DATA) {
    if (is.unsorted(-geneList))
        stop("geneList should be a decreasing sorted vector...")

    res <- GSEA_fgsea(geneList          = geneList,
                      exponent          = exponent,
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

