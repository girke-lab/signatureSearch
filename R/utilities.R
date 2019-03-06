
get_GO_Env <- function () {
  if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
  }
  get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
}

##' @importFrom S4Vectors metadata
get_organism <- function(OrgDb) {
  OrgDb <- load_OrgDb(OrgDb)
  md <- metadata(OrgDb)
  md[md[,1] == "ORGANISM", 2]
}

get_GOTERM <- function() {
    pos <- 1
    envir <- as.environment(pos)
    if (!exists(".GOTERM_Env", envir=envir)) {
        assign(".GOTERM_Env", new.env(), envir)
    }
    GOTERM_Env <- get(".GOTERM_Env", envir = envir)
    if (exists("GOTERM.df", envir = GOTERM_Env)) {
        GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
    } else {
        GOTERM.df <- toTable(GOTERM)
        assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
    }
    return(GOTERM.df)
}

get_GO2TERM_table <- function() {
    GOTERM.df <- get_GOTERM()
    GOTERM.df[, c("go_id", "Term")] %>% unique
}

add_GO_Ontology <- function(obj, GO_DATA) {
    if (is(obj, 'feaResult')) {
        obj@ontology <- 'GOALL'
    }
  
    df <- obj@result
    GO2ONT <- get("GO2ONT", envir=GO_DATA)
    df <- cbind(ONTOLOGY=GO2ONT[df$ID], df)
    obj@result <- df
    return(obj)
}

check_gene_id <- function(geneList, geneSets) {
  if (all(!names(geneList) %in% unique(unlist(geneSets)))) {
    sg <- unlist(geneSets[seq_len(10)])
    sg <- sample(sg, min(length(sg), 6))
    message("--> Expected input gene ID: ", paste0(sg, collapse=','))
    stop("--> No gene can be mapped....")
  }
}

getGeneSet <- function(USER_DATA) {
  get("PATHID2EXTID", envir = USER_DATA)
}

EXTID2TERMID <- function(gene, USER_DATA) {
  EXTID2PATHID <- get("EXTID2PATHID", envir = USER_DATA)
  
  qExtID2Path <- EXTID2PATHID[gene]
  len <- vapply(qExtID2Path, length, integer(1))
  notZero.idx <- len != 0
  qExtID2Path <- qExtID2Path[notZero.idx]
  
  return(qExtID2Path)
}

ALLEXTID <- function(USER_DATA) {
  PATHID2EXTID <- get("PATHID2EXTID", envir = USER_DATA)
  res <- unique(unlist(PATHID2EXTID))
  return(res)
}

TERMID2EXTID <- function(term, USER_DATA) {
  PATHID2EXTID <- get("PATHID2EXTID", envir = USER_DATA)
  res <- PATHID2EXTID[term]
  return(res)
}

TERM2NAME <- function(term, USER_DATA) {
  PATHID2NAME <- get("PATHID2NAME", envir = USER_DATA)
  # if (is.null(PATHID2NAME) || is.na(PATHID2NAME)) {
  #   return(as.character(term))
  # }
  return(PATHID2NAME[term])
}

get_geneSet_index <- function(geneSets, minGSSize, maxGSSize) {
  if (is.na(minGSSize) || is.null(minGSSize))
    minGSSize <- 1
  if (is.na(maxGSSize) || is.null(maxGSSize))
    maxGSSize <- Inf #.Machine$integer.max
  
  ## index of geneSets in used.
  ## logical
  geneSet_size <- vapply(geneSets, length, integer(1))
  idx <-  minGSSize <= geneSet_size & geneSet_size <= maxGSSize
  return(idx)
}

calculate_qvalue <- function(pvals) {
  if (length(pvals) == 0)
    return(numeric(0))
  
  qobj <- tryCatch(qvalue(pvals, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)
  
  if (is(qobj, "qvalue")) {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }
  return(qvalues)
}
