##' @importFrom fgsea fgsea
GSEA_fgsea <- function(geneList,
                       exponent,
                       nPerm,
                       minGSSize,
                       maxGSSize,
                       pvalueCutoff,
                       pAdjustMethod,
                       seed=FALSE,
                       USER_DATA) {
    message("preparing geneSet collections...")

    geneSets <- getGeneSet(USER_DATA)
    idx <- get_geneSet_index(geneSets, minGSSize, maxGSSize)
    geneSets <- geneSets[idx]
    
    check_gene_id(geneList, geneSets)

    message("GSEA analysis...")

    tmp_res <- fgsea(pathways=geneSets,
                 stats=geneList,
                 nperm=nPerm,
                 minSize=minGSSize,
                 maxSize=maxGSSize,
                 gseaParam=exponent,
                 nproc = 1)

    p.adj <- p.adjust(tmp_res$pval, method=pAdjustMethod)

    Description <- TERM2NAME(tmp_res$pathway, USER_DATA)

    ledge <- sapply(tmp_res$leadingEdge, paste0, collapse='/')
    ledge_rank_list <- sapply(tmp_res$leadingEdge, function(x) which(names(geneList) %in% x))
    ledge_rank <- sapply(ledge_rank_list, paste0, collapse='/')
    message("ledge_rank included")
    res <- data.frame(
        ID = as.character(tmp_res$pathway),
        Description = Description,
        setSize = tmp_res$size,
        ES = tmp_res$ES,
        NES = tmp_res$NES,
        pvalue = tmp_res$pval,
        p.adjust = p.adj,
        leadingEdge = ledge,
        ledge_rank = ledge_rank,
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
                result     = res,
                refSets   = geneSets,
                drug = names(geneList),
                universe = names(geneList)
                )
        )
    }

    row.names(res) <- res$ID
    message("done...")
    res <- as_tibble(res)
    new("feaResult",
        result     = res,
        refSets   = geneSets,
        drug = names(geneList),
        universe = names(geneList)
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
##' @param seed set seed inside the function to make result reproducible. FALSE by default.
##' @param USER_DATA annotation data
##' @return feaResult object
GSEA_internal <- function(geneList,
                 exponent,
                 nPerm,
                 minGSSize,
                 maxGSSize,
                 pvalueCutoff,
                 pAdjustMethod,
                 seed=FALSE,
                 USER_DATA) {

    if (is.unsorted(-geneList))
        geneList <- sort(geneList, decreasing = TRUE)
    
    res <- GSEA_fgsea(geneList     = geneList,
                 exponent          = exponent,
                 nPerm             = nPerm,
                 minGSSize         = minGSSize,
                 maxGSSize         = maxGSSize,
                 pvalueCutoff      = pvalueCutoff,
                 pAdjustMethod     = pAdjustMethod,
                 seed              = seed,
                 USER_DATA         = USER_DATA)
    res@organism <- "UNKNOWN"
    return(res)
}

