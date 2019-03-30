ALLEXTID <- get("ALLEXTID", envir = asNamespace("DOSE"), inherits = FALSE)
EXTID2TERMID <- get("EXTID2TERMID", envir = asNamespace("DOSE"), 
                    inherits = FALSE)
ALLEXTID <- get("ALLEXTID", envir = asNamespace("DOSE"), inherits = FALSE)
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