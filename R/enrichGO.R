##' GO Enrichment Analysis of a gene set.
##' Given a vector of genes, this function will return the enrichment GO
##' categories after FDR control.
##'
##' @param gene a vector of entrez gene id or gene SYMBOL.
##' @param OrgDb OrgDb
##' @param keytype keytype of input gene
##' @param ont One of "MF", "BP", and "CC" subontologies.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param qvalueCutoff qvalue cutoff
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param maxGSSize maximal size of genes annotated for testing
##' @param pool If ont='ALL', whether pool 3 GO sub-ontologies
##' @return A \code{feaResult} instance.
##' @seealso \code{\link{feaResult-class}}
##' @export
enrichGO <- function(gene,
                     OrgDb,
                     keytype = "ENTREZID",
                     ont="MF",
                     pvalueCutoff=0.05,
                     pAdjustMethod="BH",
                     universe,
                     qvalueCutoff = 0.2,
                     minGSSize = 10,
                     maxGSSize = 500,
                     pool=FALSE) {

    ont %<>% toupper
    ont <- match.arg(ont, c("BP", "CC", "MF", "ALL"))
    GO_DATA <- get_GO_data(OrgDb, ont, keytype)

    if (missing(universe))
        universe <- NULL

    if (ont == "ALL" && !pool) {
        lres <- lapply(c("BP", "CC", "MF"), function(ont)
            suppressMessages(enrichGO(gene, OrgDb, keytype, ont,
                     pvalueCutoff, pAdjustMethod, universe,
                     qvalueCutoff, minGSSize, maxGSSize
                     ))
            )

        lres <- lres[!sapply(lres, is.null)]
        if (length(lres) == 0)
            return(NULL)

        df <- do.call('rbind', lapply(lres, as.data.frame))
        geneSets <- lres[[1]]@geneSets
        if (length(lres) > 1) {
            for (i in 2:length(lres)) {
                geneSets <- append(geneSets, lres[[i]]@geneSets)
            }
        }
        res <- lres[[1]]
        res@result <- df
        res@refSets <- geneSets
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
    res@organism <- get_organism(OrgDb)
    res@ontology <- ont

    if (ont == "ALL") {
        res <- add_GO_Ontology(res, GO_DATA)
    }
    return(res)
}

#' get environment containing GO data
#' @title get_GO_data
#' @param OrgDb Organism annotation package, e.g., `org.Hs.eg.db` for human
#' @param ont GO ontology
#' @param keytype keytype
#' @return environment
#' @export
get_GO_data <- function(OrgDb, ont, keytype) {
  GO_Env <- get_GO_Env()
  use_cached <- FALSE
  
  if (exists("organism", envir=GO_Env, inherits=FALSE) &&
      exists("keytype", envir=GO_Env, inherits=FALSE)) {
    
    org <- get("organism", envir=GO_Env)
    kt <- get("keytype", envir=GO_Env)
    
    if (org == get_organism(OrgDb) &&
        keytype == kt &&
        exists("goAnno", envir=GO_Env, inherits=FALSE) &&
        exists("GO2TERM", envir=GO_Env, inherits=FALSE)){
      
      use_cached <- TRUE
    }
  }
  
  if (use_cached) {
    goAnno <- get("goAnno", envir=GO_Env)
  } else {
    OrgDb <- load_OrgDb(OrgDb)
    kt <- keytypes(OrgDb)
    if (! keytype %in% kt) {
      stop("keytype is not supported...")
    }
    
    kk <- keys(OrgDb, keytype=keytype)
    goAnno <- suppressMessages(
      AnnotationDbi::select(OrgDb, keys=kk, keytype=keytype,
             columns=c("GOALL", "ONTOLOGYALL")))
    
    goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ])
    
    assign("goAnno", goAnno, envir=GO_Env)
    assign("keytype", keytype, envir=GO_Env)
    assign("organism", get_organism(OrgDb), envir=GO_Env)
  }
  
  if (ont == "ALL") {
    GO2GENE <- unique(goAnno[, c(2,1)])
  } else {
    GO2GENE <- unique(goAnno[goAnno$ONTOLOGYALL == ont, c(2,1)])
  }
  
  GO_DATA <- build_Anno(GO2GENE, get_GO2TERM_table())
  
  goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
  goOnt <- goOnt.df[,2]
  names(goOnt) <- goOnt.df[,1]
  assign("GO2ONT", goOnt, envir=GO_DATA)
  return(GO_DATA)
}


get_GO_Env <- function () {
    if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
        pos <- 1
        envir <- as.environment(pos)
        assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
    }
    get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
}
