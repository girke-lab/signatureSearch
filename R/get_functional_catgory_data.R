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
#     goAnno_drug <- suppressMessages(eh[["EH3230"]])
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
    conn <- load_sqlite("EH3228")
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
    colnames(df)[1] <- keytype
    moaAnno <- df[,c("moa_id",keytype)]
    MOA_DATA <- build_Anno(moaAnno, df[,c("moa_id", "description")])
    return(MOA_DATA)
}

get_Reactome_Env <- function() {
    if (!exists(".ReactomePA_Env", envir = .GlobalEnv)) {
        assign(".ReactomePA_Env", new.env(), .GlobalEnv)
    }
    get(".ReactomePA_Env", envir= .GlobalEnv)
}

##' @importMethodsFrom AnnotationDbi as.list
##' @importFrom reactome.db reactomeEXTID2PATHID
##' @importFrom reactome.db reactomePATHID2EXTID
##' @importFrom reactome.db reactomePATHID2NAME
get_Reactome_DATA <- function(organism = "human") {
    ReactomePA_Env <- get_Reactome_Env()
    
    if (exists("organism", envir=ReactomePA_Env, inherits = FALSE)) {
        org <- get("organism", envir=ReactomePA_Env)
        if (org == organism &&
            exists("PATHID2EXTID", envir = ReactomePA_Env) &&
            exists("EXTID2PATHID", envir = ReactomePA_Env) &&
            exists("PATHID2NAME",  envir = ReactomePA_Env)) {
            
            ## use_cached
            return(ReactomePA_Env)
        }
    }
    
    ALLEG <- getALLEG(organism)
    
    EXTID2PATHID <- as.list(reactomeEXTID2PATHID)
    EXTID2PATHID <- EXTID2PATHID[names(EXTID2PATHID) %in% ALLEG]
    
    PATHID2EXTID <- as.list(reactomePATHID2EXTID) ## also contains reactions
    
    PATHID2NAME <- as.list(reactomePATHID2NAME)
    PI <- names(PATHID2NAME)
    PATHID2NAME <- lapply(PATHID2NAME, function(x) x[1])
    names(PATHID2NAME) <- PI
    
    PATHID2EXTID <- PATHID2EXTID[names(PATHID2EXTID) %in% names(PATHID2NAME)]
    PATHID2EXTID <- PATHID2EXTID[names(PATHID2EXTID) %in% unique(unlist(EXTID2PATHID))]
    PATHID2EXTID <- lapply(PATHID2EXTID, function(x) intersect(x, ALLEG))
    
    PATHID2NAME <- PATHID2NAME[names(PATHID2NAME) %in% names(PATHID2EXTID)]
    
    PATHID2NAME <- unlist(PATHID2NAME)
    PATHID2NAME <- gsub("^\\w+\\s\\w+:\\s+", "", PATHID2NAME) # remove leading spaces
    
    assign("PATHID2EXTID", PATHID2EXTID, envir=ReactomePA_Env)
    assign("EXTID2PATHID", EXTID2PATHID, envir=ReactomePA_Env)
    assign("PATHID2NAME", PATHID2NAME, envir=ReactomePA_Env)
    return(ReactomePA_Env)
}

##' get all entrez gene ID of a specific organism
##'
##' @title getALLEG
##' @param organism one of "human", "rat", "mouse", "celegans", "yeast", 
##' "zebrafish", "fly".
##' @return entrez gene ID vector
##' @importMethodsFrom AnnotationDbi keys
##' @author Yu Guangchuang
getALLEG <- function(organism) {
    annoDb <- getDb(organism)
    annoDb <- load_OrgDb(annoDb)
    eg <- keys(annoDb, keytype="ENTREZID")
    return(eg)
}

##' mapping organism name to annotationDb package name
##'
##' @title getDb
##' @param organism one of supported organism
##' @return annotationDb name
##' @author Yu Guangchuang
getDb <- function(organism) {
    annoDb <- switch(organism,
                     anopheles   = "org.Ag.eg.db",
                     arabidopsis = "org.At.tair.db",
                     bovine      = "org.Bt.eg.db",
                     canine      = "org.Cf.eg.db",
                     celegans    = "org.Ce.eg.db",
                     chicken     = "org.Gg.eg.db",
                     chimp       = "org.Pt.eg.db",
                     coelicolor  = "org.Sco.eg.db", 
                     ecolik12    = "org.EcK12.eg.db",
                     ecsakai     = "org.EcSakai.eg.db",
                     fly         = "org.Dm.eg.db",
                     gondii      = "org.Tgondii.eg.db",
                     human       = "org.Hs.eg.db",
                     malaria     = "org.Pf.plasmo.db",
                     mouse       = "org.Mm.eg.db",
                     pig         = "org.Ss.eg.db",
                     rat         = "org.Rn.eg.db",
                     rhesus      = "org.Mmu.eg.db",
                     xenopus     = "org.Xl.eg.db",
                     yeast       = "org.Sc.sgd.db",
                     zebrafish   = "org.Dr.eg.db",
    )
    return(annoDb)
}

