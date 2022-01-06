ALLEXTID <- get("ALLEXTID", envir = asNamespace("DOSE"), inherits = FALSE)
EXTID2TERMID <- get("EXTID2TERMID", envir = asNamespace("DOSE"), 
                    inherits = FALSE)
EXTID2NAME <- get("EXTID2NAME", envir = asNamespace("DOSE"), inherits = FALSE)
TERM2NAME <- get("TERM2NAME", envir = asNamespace("DOSE"), inherits = FALSE)
TERMID2EXTID <- get("TERMID2EXTID", envir = asNamespace("DOSE"), 
                    inherits = FALSE)
build_Anno <- get("build_Anno", envir = asNamespace("DOSE"), inherits = FALSE)
calculate_qvalue <- get("calculate_qvalue", envir = asNamespace("DOSE"), 
                        inherits = FALSE)
geneSet_filter <- get("geneSet_filter", envir = asNamespace("DOSE"), 
                      inherits = FALSE)
get_geneSet_index <- get("get_geneSet_index", envir = asNamespace("DOSE"), 
                         inherits = FALSE)
get_organism <- get("get_organism", envir = asNamespace("DOSE"), 
                    inherits = FALSE)

add_GO_Ontology <- get("add_GO_Ontology", 
                       envir = asNamespace("clusterProfiler"), inherits = FALSE)
get_GO2TERM_table <- get("get_GO2TERM_table", 
                       envir = asNamespace("clusterProfiler"), inherits = FALSE)
get_GO_Env <- get("get_GO_Env", 
                       envir = asNamespace("clusterProfiler"), inherits = FALSE)
organismMapper <- get("organismMapper", 
                       envir = asNamespace("clusterProfiler"), inherits = FALSE)
prepare_KEGG <- get("prepare_KEGG", 
                       envir = asNamespace("clusterProfiler"), inherits = FALSE)

validh5 <- function(ehid){
    h5path <- eh[[ehid]]
    tryCatch(h5ls(h5path), error=function(e){
        unlink(h5path)
        h5path <- eh[[ehid]]
    })
    return(h5path)
}

determine_refdb <- function(refdb){
    if(refdb=="cmap") return(validh5("EH3223"))
    if(refdb=="cmap_expr") return(validh5("EH3224"))
    if(refdb=="lincs") return(validh5("EH3226"))
    if(refdb=="lincs_expr") return(validh5("EH3227"))
    if(refdb=="lincs2") return(validh5("EH7297"))
    return(refdb)
}

load_sqlite <- function(ehid){
    path <- suppressMessages(eh[[ehid]])
    conn <- tryCatch(dbConnect(SQLite(), path), error=function(e){
        unlink(path)
        path <- eh[[ehid]]
        dbConnect(SQLite(), path)
    })
    return(conn)
}