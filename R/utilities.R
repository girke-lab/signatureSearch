
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

# ##' convert goid to descriptive term
# ##'
# ##'
# ##' @title go2term
# ##' @param goid a vector of GO IDs
# ##' @return data.frame
# ##' @export
# ##' @author Guangchuang Yu
# go2term <- function(goid) {
#     GO2TERM <- get_GO2TERM_table()
#     res <- GO2TERM[GO2TERM[,1] %in% goid, ]
#     rownames(res) <- NULL
#     return(res)
# }
# 
# ##' convert goid to ontology (BP, CC, MF)
# ##'
# ##'
# ##' @title go2ont
# ##' @param goid a vector of GO IDs
# ##' @return data.frame
# ##' @export
# ##' @author Guangchuang Yu
# go2ont <- function(goid) {
#     GO2Ontology <- get_GO2Ontology_table()
#     res <- GO2Ontology[GO2Ontology[,1] %in% goid,]
#     rownames(res) <- NULL
#     return(res)
# }


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

# get_GO2Ontology_table <- function() {
#     GOTERM.df <- get_GOTERM()
#     GOTERM.df[, c("go_id", "Ontology")] %>% unique
# }
# 
# 
# excludeGOlevel <- function(x, ont, level) {
#     lv <- unlist(lapply(level, getGOLevel, ont=ont))
#     x <- excludeGOterm(x, lv)
#     return(x)
# }
# 
# excludeGOterm <- function(x, term) {
#     if ( is(x, "enrichResult") ) {
#         x@result <- x@result[! x@result[, "ID"] %in% term, ]
#     } else if ( is(x, "compareClusterResult") ) {
#         x@compareClusterResult <- x@compareClusterResult[! x@compareClusterResult[, "ID"] %in% term, ]
#     } else {
#         stop("x should be one of enrichResult of compareClusterResult...")
#     }
#     return(x)
# }
# 
# keepGOlevel <- function(x, ont, level) {
#     lv <- unlist(lapply(level, getGOLevel, ont=ont))
#     x <- keepGOterm(x, lv)
#     return(x)
# }
# 
# keepGOterm <- function(x, term) {
#     if ( is(x, "enrichResult") ) {
#         x@result <- x@result[x@result[, "ID"] %in% term, ]
#     } else if ( is(x, "compareClusterResult") ) {
#         x@compareClusterResult <- x@compareClusterResult[x@compareClusterResult[, "ID"] %in% term, ]
#     } else {
#         stop("x should be one of enrichResult of compareClusterResult...")
#     }
#     return(x)
# }
# 
# ##' query GOIDs at a specific level.
# ##'
# ##'
# ##' @title get GOIDs at a specific level
# ##' @param ont Ontology
# ##' @param level GO level
# ##' @return a vector of GOIDs
# ##' @importFrom GO.db GOBPCHILDREN
# ##' @importFrom GO.db GOCCCHILDREN
# ##' @importFrom GO.db GOMFCHILDREN
# ##' @importMethodsFrom AnnotationDbi mget
# ##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
# getGOLevel <- function(ont, level) {
#     switch(ont,
#            MF = {
#                topNode <- "GO:0003674"
#                Children <- GOMFCHILDREN
#            },
#            BP = {
#                topNode <- "GO:0008150"
#                Children <- GOBPCHILDREN
#            },
#            CC = {
#                topNode <- "GO:0005575"
#                Children <- GOCCCHILDREN
#            }
#            )
# 
#     max_level <- max(level)
#     if (any(level == 1)) {
#         all_nodes <- topNode
#     } else {
#         all_nodes <- c()
#     }
# 
#     Node <- topNode
#     for (i in seq_len(max_level-1)) {
#         Node <- mget(Node, Children, ifnotfound=NA)
#         Node <- unique(unlist(Node))
#         Node <- as.vector(Node)
#         Node <- Node[!is.na(Node)]
#         if ((i+1) %in% level) {
#             all_nodes <- c(all_nodes, Node)
#         }
#     }
#     return(all_nodes)
# }
# 
# 

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
    sg <- unlist(geneSets[1:10])
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
  len <- sapply(qExtID2Path, length)
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
  if (is.null(PATHID2NAME) || is.na(PATHID2NAME)) {
    return(as.character(term))
  }
  return(PATHID2NAME[term])
}

get_geneSet_index <- function(geneSets, minGSSize, maxGSSize) {
  if (is.na(minGSSize) || is.null(minGSSize))
    minGSSize <- 1
  if (is.na(maxGSSize) || is.null(maxGSSize))
    maxGSSize <- Inf #.Machine$integer.max
  
  ## index of geneSets in used.
  ## logical
  geneSet_size <- sapply(geneSets, length)
  idx <-  minGSSize <= geneSet_size & geneSet_size <= maxGSSize
  return(idx)
}

# 
# get_go_ontology <- function(x) {
#     if (is(x, "compareClusterResult")) {
#         if (x@fun != "enrichGO" && x@fun != "groupGO") {
#             stop("simplify only work for GO...")
#         }
#         ont <- x@.call$ont
#         if (is.null(ont) || class(ont) != "character") {
#             ## should be "MF", default value of enrichGO
#             ## it's safe to determine from the output
#             ont <- x@compareClusterResult$ID[1] %>% GOTERM[[.]] %>% Ontology
#         }
#     } else if (is(x, "enrichResult")) {
#         if (!x@ontology %in% c("BP", "MF", "CC"))
#             stop("ontology should be one of 'MF', 'BP', 'CC'...")
# 
#         ont <- x@ontology
#     } else {
#         stop("x should be an instance of 'enrichResult' or 'compareClusterResult'...")
#     }
# 
#     return(ont)
# }
