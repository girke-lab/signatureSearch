##' @importFrom AnnotationDbi keys
##' @importFrom AnnotationDbi select
##' @importFrom AnnotationDbi keytypes
##' @importFrom AnnotationDbi toTable
##' @importFrom dplyr left_join
##' @importFrom dplyr as_tibble
##' @importFrom dplyr filter
##' @importFrom dplyr distinct
##' @importFrom magrittr %>%
##' @importFrom RSQLite dbConnect
##' @importFrom RSQLite dbGetQuery
##' @importFrom RSQLite SQLite
##' @importFrom RSQLite dbDisconnect
##' @importFrom utils download.file
##' @importFrom stats na.omit
##' @importFrom GOSemSim load_OrgDb
# get_GO_data_drug <- function(OrgDb, ont, keytype) {
#     GO_Env <- get_GO_Env()
#     use_cached <- FALSE
#     
#     if (exists("organism", envir=GO_Env, inherits=FALSE) &&
#         exists("keytype", envir=GO_Env, inherits=FALSE)) {
#         
#         org <- get("organism", envir=GO_Env)
#         kt <- get("keytype", envir=GO_Env)
#         
#         if (org == get_organism(OrgDb) &&
#             keytype == kt &&
#             exists("goAnno", envir=GO_Env, inherits=FALSE) &&
#             exists("GO2TERM", envir=GO_Env, inherits=FALSE)){
#             
#             use_cached <- TRUE
#         }
#     }
#     
#     if (use_cached) {
#         goAnno <- get("goAnno", envir=GO_Env)
#     } else {
#         OrgDb <- load_OrgDb(OrgDb)
#         kt <- keytypes(OrgDb)
#         if (! keytype %in% kt) {
#             stop("keytype is not supported...")
#         }
#         
#         kk <- keys(OrgDb, keytype=keytype)
#         goAnno <- suppressMessages(
#             AnnotationDbi::select(OrgDb, keys=kk, keytype=keytype,
#                                   columns=c("GOALL", "ONTOLOGYALL")))
#         
#         goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ])
#         
#         assign("goAnno", goAnno, envir=GO_Env)
#         assign("keytype", keytype, envir=GO_Env)
#         assign("organism", get_organism(OrgDb), envir=GO_Env)
#     }
#     
#     # download goAnno_drug.rds 
#     goAnno_drug <- suppressMessages(ah[["AH69085"]])
#     ## "drug_name" in goAnno_drug are all lowercase
#     
#     if (ont == "ALL") {
#         GO2GENE <- goAnno_drug[,c("GOALL","drug_name")]
#     } else {
#         GO2GENE <- goAnno_drug[goAnno_drug$ONTOLOGYALL == ont, 
#                                c("GOALL","drug_name")]
#     }
#     
#     GO_DATA <- build_Anno(GO2GENE, get_GO2TERM_table())
#     
#     goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
#     goOnt <- goOnt.df[,2]
#     names(goOnt) <- goOnt.df[,1]
#     assign("GO2ONT", goOnt, envir=GO_DATA)
#     return(GO_DATA)
# }

##' @importFrom dplyr left_join
##' @importFrom dplyr as_tibble
##' @importFrom dplyr filter
##' @importFrom dplyr distinct
##' @importFrom magrittr %>%
##' @importFrom RSQLite dbConnect
##' @importFrom RSQLite dbGetQuery
##' @importFrom RSQLite SQLite
##' @importFrom RSQLite dbDisconnect
##' @importFrom utils download.file
##' @importFrom clusterProfiler download_KEGG
prepare_KEGG_drug <- function(species, KEGG_Type="KEGG", keyType="kegg") {
    kegg <- clusterProfiler::download_KEGG(species, KEGG_Type, keyType)
    # get dtlink_entrez
    conn <- load_sqlite("AH69083")
    dtlink_entrez <- dbGetQuery(conn, 'SELECT * FROM dtlink_entrez')
    dbDisconnect(conn)
    
    # map entrez id to drug name
    keggpath2entrez <- kegg$KEGGPATHID2EXTID
    keggpath2entrez <- left_join(as_tibble(keggpath2entrez), 
                                 as_tibble(dtlink_entrez), 
                                 by = c("to"="ENTREZID")) %>% 
        filter(!is.na(drug_name)) %>% 
        distinct(from, drug_name, .keep_all = TRUE)
    keggpath2entrez <- as.data.frame(keggpath2entrez)[,c("from","drug_name")]
    kegg$KEGGPATHID2EXTID <- keggpath2entrez
    build_Anno(kegg$KEGGPATHID2EXTID, kegg$KEGGPATHID2NAME)
}

## Build MOA_DATA environment
get_MOA_data <- function(moa_list, keytype="drug_name") {
    moa_list2 <- moa_list
    names(moa_list2) <- paste0("MOA:", sprintf("%04d", seq_along(moa_list2)))
    df <- data.frame(unname=unlist(moa_list), 
                     moa_id=rep(names(moa_list2), lengths(moa_list2)),
                     description=rep(names(moa_list), lengths(moa_list)), 
                     stringsAsFactors=FALSE, row.names=NULL)
    colnames(df)[1] = keytype
    moaAnno <- df[,c("moa_id",keytype)]
    MOA_DATA <- build_Anno(moaAnno, df[,c("moa_id", "description")])
    return(MOA_DATA)
}

## get rid of "Undefined global functions or variables" note
## cell = drug_name = from = . = NULL
globalVariables(c("cell", "drug_name", "from")) 
