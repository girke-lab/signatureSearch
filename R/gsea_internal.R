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

    geneSets <- getGeneSet(USER_DATA)
    check_gene_id(geneList, geneSets)

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
    qvalues <- calculate_qvalue(tmp_res$pval)

    Description <- TERM2NAME(tmp_res$pathway, USER_DATA)

    params <- list(pvalueCutoff = pvalueCutoff,
                   nPerm = nPerm,
                   pAdjustMethod = pAdjustMethod,
                   exponent = exponent,
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize
                   )
    ledge <- vapply(tmp_res$leadingEdge, paste0, collapse='/', FUN.VALUE = character(1))
    ledge_rank <- lapply(tmp_res$leadingEdge, function(x) match(x, names(geneList)))
    ledge_rank2 <- vapply(ledge_rank, paste, collapse="/", FUN.VALUE = character(1))
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
        message("no term enriched under specific pvalueCutoff...")
        return(
          new("feaResult",
              result    = as_tibble(res),
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

##' generic function for gene set enrichment analysis
##'
##'
##' @title GSEA_internal
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
##' @param by one of 'fgsea' or 'DOSE'
##' @return gseaResult object
##' @author Yu Guangchuang
GSEA_internal <- function(geneList,
                 exponent,
                 nPerm,
                 minGSSize,
                 maxGSSize,
                 pvalueCutoff,
                 pAdjustMethod,
                 verbose=TRUE,
                 seed=FALSE,
                 USER_DATA,
                 by="fgsea") {

    by <- match.arg(by, c("fgsea", "DOSE"))
    if (is.unsorted(-geneList))
        stop("geneList should be a decreasing sorted vector...")
    if (by == 'fgsea') {
        .GSEA <- GSEA_fgsea
    } else {
        .GSEA <- GSEA_fgsea
    }

    res <- .GSEA(geneList          = geneList,
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

