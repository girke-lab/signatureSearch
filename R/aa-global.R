
ALLEXTID <- function (USER_DATA){
  if (inherits(USER_DATA, "environment")) {
    PATHID2EXTID <- get("PATHID2EXTID", envir = USER_DATA)
    res <- unique(unlist(PATHID2EXTID))
  }
  else if (inherits(USER_DATA, "GSON")) {
    gsid2gene <- USER_DATA@gsid2gene
    res <- unique(gsid2gene$gene)
  }
  else {
    stop("not supported")
  }
  return(res)
}

EXTID2TERMID <- function (gene, USER_DATA){
  if (inherits(USER_DATA, "environment")) {
    EXTID2PATHID <- get("EXTID2PATHID", envir = USER_DATA)
    qExtID2Path <- EXTID2PATHID[gene]
  }
  else if (inherits(USER_DATA, "GSON")) {
    gsid2gene <- USER_DATA@gsid2gene
    qExtID2Path <- setNames(lapply(gene, function(x) {
      subset(gsid2gene, gsid2gene$gene == x)[["gsid"]]
    }), gene)
  }
  else {
    stop("not supported")
  }
  len <- sapply(qExtID2Path, length)
  notZero.idx <- len != 0
  qExtID2Path <- qExtID2Path[notZero.idx]
  return(qExtID2Path)
}

#' @importFrom AnnotationDbi select
EXTID2NAME <- function(OrgDb, geneID, keytype){
  OrgDb <- load_OrgDb(OrgDb)
  kt <- keytypes(OrgDb)
  if (!keytype %in% kt) {
    stop("keytype is not supported...")
  }
  gn.df <- suppressMessages(AnnotationDbi::select(OrgDb, keys = geneID, keytype = keytype, 
                                   columns = "SYMBOL"))
  gn.df <- unique(gn.df)
  colnames(gn.df) <- c("GeneID", "SYMBOL")
  unmap_geneID <- geneID[!geneID %in% gn.df$GeneID]
  if (length(unmap_geneID) != 0) {
    unmap_geneID.df <- data.frame(GeneID = unmap_geneID, SYMBOL = unmap_geneID)
    gn.df <- rbind(gn.df, unmap_geneID.df)
  }
  gn <- gn.df$SYMBOL
  names(gn) <- gn.df$GeneID
  return(gn)
}

TERM2NAME <- function (term, USER_DATA) {
  if (inherits(USER_DATA, "environment")) {
    PATHID2NAME <- get("PATHID2NAME", envir = USER_DATA)
    if (is.null(PATHID2NAME) || all(is.na(PATHID2NAME))) {
      return(as.character(term))
    }
    return(PATHID2NAME[term])
  }
  else if (inherits(USER_DATA, "GSON")) {
    gsid2name <- USER_DATA@gsid2name
    res <- setNames(vapply(term, function(x) {
      subset(gsid2name, gsid2name$gsid == x)[["name"]]
    }, character(1)), term)
    return(res)
  }
  return(as.character(term))
}

TERMID2EXTID <- function (term, USER_DATA){
  if (inherits(USER_DATA, "environment")) {
    PATHID2EXTID <- get("PATHID2EXTID", envir = USER_DATA)
    res <- PATHID2EXTID[term]
  }
  else if (inherits(USER_DATA, "GSON")) {
    gsid2gene <- USER_DATA@gsid2gene
    res <- setNames(lapply(term, function(x) {
      subset(gsid2gene, gsid2gene$gsid == x)[["gene"]]
    }), term)
  }
  else {
    stop("not supported")
  }
  return(res)
}

build_Anno <- function (path2gene, path2name){
  if (!exists(".Anno_clusterProfiler_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".Anno_clusterProfiler_Env", new.env(), envir = envir)
  }
  Anno_clusterProfiler_Env <- get(".Anno_clusterProfiler_Env", 
                                  envir = .GlobalEnv)
  if (inherits(path2gene[[2]], "list")) {
    path2gene <- cbind(rep(path2gene[[1]], times = vapply(path2gene[[2]], 
                                                          length, numeric(1))), unlist(path2gene[[2]]))
  }
  path2gene <- as.data.frame(path2gene)
  path2gene <- path2gene[!is.na(path2gene[, 1]), ]
  path2gene <- path2gene[!is.na(path2gene[, 2]), ]
  path2gene <- unique(path2gene)
  PATHID2EXTID <- split(as.character(path2gene[, 2]), as.character(path2gene[, 
                                                                             1]))
  EXTID2PATHID <- split(as.character(path2gene[, 1]), as.character(path2gene[, 
                                                                             2]))
  assign("PATHID2EXTID", PATHID2EXTID, envir = Anno_clusterProfiler_Env)
  assign("EXTID2PATHID", EXTID2PATHID, envir = Anno_clusterProfiler_Env)
  if (missing(path2name) || is.null(path2name) || all(is.na(path2name))) {
    assign("PATHID2NAME", NULL, envir = Anno_clusterProfiler_Env)
  }
  else {
    path2name <- as.data.frame(path2name)
    path2name <- path2name[!is.na(path2name[, 1]), ]
    path2name <- path2name[!is.na(path2name[, 2]), ]
    path2name <- unique(path2name)
    PATH2NAME <- as.character(path2name[, 2])
    names(PATH2NAME) <- as.character(path2name[, 1])
    assign("PATHID2NAME", PATH2NAME, envir = Anno_clusterProfiler_Env)
  }
  return(Anno_clusterProfiler_Env)
}

calculate_qvalue <- function(pvals){
  if (length(pvals) == 0) 
    return(numeric(0))
  qobj <- tryCatch(qvalue(pvals, lambda = 0.05, pi0.method = "bootstrap"), 
                   error = function(e) NULL)
  if (inherits(qobj, "qvalue")) {
    qvalues <- qobj$qvalues
  }
  else {
    qvalues <- NA
  }
  return(qvalues)
}

geneSet_filter <- function (geneSets, geneList, minGSSize, maxGSSize){
  geneSets <- sapply(geneSets, intersect, names(geneList))
  gs.idx <- get_geneSet_index(geneSets, minGSSize, maxGSSize)
  nGeneSet <- sum(gs.idx)
  if (nGeneSet == 0) {
    msg <- paste0("No gene set have size between [", minGSSize, 
                  ", ", maxGSSize, "]...")
    message(msg)
    message("--> return NULL...")
    return(NULL)
  }
  geneSets[gs.idx]
}

get_geneSet_index <- function (geneSets, minGSSize, maxGSSize){
  if (is.na(minGSSize) || is.null(minGSSize)) 
    minGSSize <- 1
  if (is.na(maxGSSize) || is.null(maxGSSize)) 
    maxGSSize <- Inf
  geneSet_size <- sapply(geneSets, length)
  idx <- minGSSize <= geneSet_size & geneSet_size <= maxGSSize
  return(idx)
}

#' @importFrom AnnotationDbi species
get_organism <- function(OrgDb){
  OrgDb <- load_OrgDb(OrgDb)
  AnnotationDbi::species(OrgDb)
}

add_GO_Ontology <- function (obj, GO_DATA){
  if (is(obj, "gseaResult")) {
    obj@setType <- "GOALL"
  }
  else if (is(obj, "enrichResult")) {
    obj@ontology <- "GOALL"
  }
  df <- obj@result
  GO2ONT <- get("GO2ONT", envir = GO_DATA)
  df <- cbind(ONTOLOGY = GO2ONT[df$ID], df)
  obj@result <- df
  return(obj)
}

#' @import annotate 
get_GO2TERM_table <- function(){
  GOTERM.df <- get_GOTERM()
  GOTERM.df[, c("go_id", "Term")] %>% unique
}

get_GO_Env <- function(){
  if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".GO_clusterProfiler_Env", new.env(), envir = envir)
  }
  get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
}

organismMapper <- function(organism){
  if (organism == "anopheles") {
    species <- "aga"
  } else if (organism == "arabidopsis") {
    species <- "ath"
  } else if (organism == "bovine") {
    species <- "bta"
  } else if (organism == "canine") {
    species <- "cfa"
  } else if (organism == "chicken") {
    species <- "gga"
  } else if (organism == "chipm") {
    species <- "ptr"
  } else if (organism == "ecolik12") {
    species <- "eco"
  } else if (organism == "ecsakai") {
    species <- "ecs"
  } else if (organism == "fly") {
    species <- "dme"
  } else if (organism == "human") {
    species <- "hsa"
  } else if (organism == "malaria") {
    species <- "pfa"
  } else if (organism == "mouse") {
    species <- "mmu"
  } else if (organism == "pig") {
    species <- "ssc"
  } else if (organism == "rat") {
    species <- "rno"
  } else if (organism == "rhesus") {
    species <- "mcc"
  } else if (organism == "worm" || organism == "celegans") {
    species <- "cel"
  } else if (organism == "xenopus") {
    species <- "xla"
  } else if (organism == "yeast") {
    species <- "sce"
  } else if (organism == "zebrafish") {
    species <- "dre"
  } else {
    species <- organism
  }
  return(species)
}

#' @import ExperimentHub
validh5 <- function(ehid){
    eh <- suppressMessages(ExperimentHub())
    h5path <- eh[[ehid]]
    tryCatch(h5ls(h5path), error=function(e){
        unlink(h5path)
        h5path <- eh[[ehid]]
    })
    return(h5path)
}

determine_refdb <- function(refdb){
    eh <- suppressMessages(ExperimentHub())
    if(refdb=="cmap") return(validh5("EH3223"))
    if(refdb=="cmap_expr") return(validh5("EH3224"))
    if(refdb=="lincs") return(validh5("EH3226"))
    if(refdb=="lincs_expr") return(validh5("EH3227"))
    if(refdb=="lincs2") return(validh5("EH7297"))
    return(refdb)
}

load_sqlite <- function(ehid){
    eh <- suppressMessages(ExperimentHub())
    path <- suppressMessages(eh[[ehid]])
    conn <- tryCatch(dbConnect(SQLite(), path), error=function(e){
        unlink(path)
        path <- eh[[ehid]]
        dbConnect(SQLite(), path)
    })
    return(conn)
}

#' @importFrom BiocGenerics fileName
validLoad <- function(ehid){
    eh <- suppressMessages(ExperimentHub())
    tryCatch(suppressMessages(eh[[ehid]]), 
             error=function(e){
                 unlink(fileName(eh[ehid]))
                 eh[[ehid]]})
}

# GO_DATA <- get_GO_data(OrgDb, ont, keytype="SYMBOL")
# download GO_DATA.rds from AnnotationHub to save time by avoiding 
# building GO_DATA from scratch
GO_DATA <- validLoad("EH3231")

# GO_DATA_drug <- get_GO_data_drug(OrgDb = "org.Hs.eg.db", 
#                                  ont, keytype="SYMBOL")
# download GO_DATA_drug.rds 
GO_DATA_drug <- validLoad("EH3232")
