##' @importFrom qvalue qvalue
##' @importFrom stats phyper
##' @importFrom stats p.adjust

## Interal method for enrichment analysis by using the hypergeometric model, 
## which also support query gene set with duplications
enricher_internal <- function(gene,
                              pvalueCutoff,
                              pAdjustMethod="BH",
                              universe = NULL,
                              minGSSize=10,
                              maxGSSize=500,
                              qvalueCutoff=0.2,
                              USER_DATA){

  ## query external ID to Term ID
  gene <- as.character(gene)
  qExtID2TermID <- DOSE:::EXTID2TERMID(gene, USER_DATA)
  qTermID <- unlist(qExtID2TermID)
  if (is.null(qTermID)) {
      message("--> No gene can be mapped....")

      p2e <- get("PATHID2EXTID", envir=USER_DATA)
      sg <- unlist(p2e[seq_len(10)])
      sg <- sample(sg, min(length(sg), 6))
      message("--> Expected input gene ID: ", paste0(sg, collapse=','))

      message("--> return NULL...")
      return(NULL)
  }

  ## Term ID -- query external ID association list.
  qExtID2TermID.df <- data.frame(extID=rep(names(qExtID2TermID),
                                           times=lapply(qExtID2TermID, length)),
                                 termID=qTermID)
  qExtID2TermID.df <- unique(qExtID2TermID.df)

  qTermID2ExtID <- split(as.character(qExtID2TermID.df$extID), 
                         as.character(qExtID2TermID.df$termID))

  ## Get all the genes that have GO annotation, intersect with universe, 
  ## get extID as universe
  extID <- DOSE:::ALLEXTID(USER_DATA)
  if (missing(universe))
      universe <- NULL
  if(!is.null(universe)) {
      extID <- intersect(extID, universe)
  }
  
  ## get intersect of query genes and universe, as new query genes
  qGene <- gene[gene %in% extID]

  ## Term ID annotate query external ID
  qTermID <- unique(names(qTermID2ExtID))


  termID2ExtID <- DOSE:::TERMID2EXTID(qTermID, USER_DATA)
  termID2ExtID <- lapply(termID2ExtID, intersect, extID)

  idx <- DOSE:::get_geneSet_index(termID2ExtID, minGSSize, maxGSSize)

  if (sum(idx) == 0) {
      msg <- paste("No gene set have size >", minGSSize, "...")
      message(msg)
      message("--> return NULL...")
      return (NULL)
  }

  termID2ExtID <- termID2ExtID[idx]
  geneSets <- termID2ExtID
  qTermID2ExtID <- qTermID2ExtID[idx]
  qTermID <- unique(names(qTermID2ExtID))

  ## prepare parameter for hypergeometric test
  
  # N: all the balls in the urn.
  N <- length(extID)
  
  # k: White balls drawn. number of overlapped genes in qGene and genes in 
  # every GO term of termID2ExtID
  k <- vapply(termID2ExtID, function(x) sum(qGene %in% x), 
              FUN.VALUE = integer(1))

  # n: balls drawn. length of qGene
  n <- length(qGene)
  
  # M:  White balls.
  M <- vapply(termID2ExtID, length, FUN.VALUE = integer(1))

  args.df <- data.frame(numWdrawn=k-1, ## White balls drawn
                        numW=M,        ## White balls
                        numB=N-M,      ## Black balls
                        numDrawn=n)    ## balls drawn


  ## calcute pvalues based on hypergeometric model
  pvalues <- apply(args.df, 1, function(n)
                   phyper(n[1], n[2], n[3], n[4], lower.tail=FALSE)
                   )

  ## gene ratio and background ratio
  GeneRatio <- apply(data.frame(a=k, b=n), 1, function(x)
                     paste(x[1], "/", x[2], sep="", collapse="")
                     )
  BgRatio <- apply(data.frame(a=M, b=N), 1, function(x)
                   paste(x[1], "/", x[2], sep="", collapse="")
                   )

  Over <- data.frame(ID = as.character(qTermID),
                     GeneRatio = GeneRatio,
                     BgRatio = BgRatio,
                     pvalue = pvalues,
                     stringsAsFactors = FALSE)

  p.adj <- p.adjust(Over$pvalue, method=pAdjustMethod)
  qobj <- tryCatch(qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap"), 
                   error=function(e) NULL)

  if (is(qobj, "qvalue")) {
      qvalues <- qobj$qvalues
  } else {
      qvalues <- NA
  }
  
  ## Obtain genes matching at GO nodes of termID2ExtID
  geneID <- lapply(termID2ExtID, function(x) { qGene[qGene %in% x] } )
  geneID <- vapply(geneID, function(x) { paste(x, collapse="/") }, 
                   FUN.VALUE = character(1) )
  geneID[geneID==""] <- "NA"
  geneID <- geneID[qTermID]
  Over <- data.frame(Over,
                     p.adjust = p.adj,
                     qvalue = qvalues,
                     geneID = geneID,
                     Count = k,
                     stringsAsFactors = FALSE)

  Description <- DOSE:::TERM2NAME(qTermID, USER_DATA)

  if (length(qTermID) != length(Description)) {
      idx <- qTermID %in% names(Description)
      Over <- Over[idx,]
  }
  Over$Description <- Description
  nc <- ncol(Over)
  Over <- Over[, c(1,nc, 2:(nc-1))]


  Over <- Over[order(pvalues),]

  Over <- Over[ Over$pvalue <= pvalueCutoff, ]
  Over <- Over[ Over$p.adjust <= pvalueCutoff, ]
  if (! any(is.na(Over$qvalue))) {
      Over <- Over[ Over$qvalue <= qvalueCutoff, ]
  }

  Over$ID <- as.character(Over$ID)
  Over$Description <- as.character(Over$Description)

  row.names(Over) <- as.character(Over$ID)
  Over <- as_tibble(Over)
  x <- new("feaResult",
           result         = Over,
           drugs          = as.character(gene),
           targets        = as.character(gene),
           #universe       = extID,
           #refSets        = geneSets,
           organism       = "UNKNOWN",
           ontology       = "UNKNOWN"
           )
  return (x)
}
