##' Drug GO Enrichment Analysis of drug sets.
##' Given a vector of drugs, this function will return the enriched GO categories with hypergeometric test after FDR control.
##'
##'
##' @param drugs query drug set used to do drug set enrichment analysis (DSEA), Can be top ranking drugs in GESS result. 
##' @param type one of "GO" or "KEGG"
##' @param ont One of "MF", "BP", and "CC" or "ALL".
##' @param pvalueCutoff Cutoff value of p value.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param qvalueCutoff q value cutoff
##' @param minGSSize minimal size of drugs annotated by Ontology term for testing after drug-to-GO category mapping.
##' @param maxGSSize maximal size of drugs annotated for testing
##' @return \code{feaResult} object
##' @import org.Hs.eg.db
##' @importFrom magrittr %<>%
##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importMethodsFrom AnnotationDbi mget
##' @seealso \code{\link{feaResult-class}}
##' @export
##' @author Yuzhu Duan
dsea_hyperG <- function(drugs,
                     type = "GO",
                     ont="BP",
                     pvalueCutoff=0.05,
                     pAdjustMethod="BH",
                     qvalueCutoff = 0.2,
                     minGSSize = 10,
                     maxGSSize = 500) {
    drugs <- as.character(unique(tolower(drugs)))
    if(type == "GO"){
      ont %<>% toupper
      ont <- match.arg(ont, c("BP", "CC", "MF", "ALL"))

      # GO_DATA_drug <- get_GO_data_drug(OrgDb = "org.Hs.eg.db", ont, keytype="SYMBOL")
      # download GO_DATA_drug and save it locally to save time
      ext_path <- system.file("extdata", package="signatureSearch")
      godata_drug_path <- paste0(ext_path,"/GO_DATA_drug.rds")
      if(file.exists(godata_drug_path)){
        GO_DATA_drug <- readRDS(godata_drug_path)
      } else {
        download.file("http://biocluster.ucr.edu/~yduan004/fea/GO_DATA_drug.rds", godata_drug_path, quiet = TRUE)
        GO_DATA_drug <- readRDS(godata_drug_path)
      }
      
      # get all the drugs in the corresponding annotation system as universe
      ext2path <- get("EXTID2PATHID", envir = GO_DATA_drug)
      universe = names(ext2path)
      
      res <- enricher_internal(gene=drugs,
                               pvalueCutoff=pvalueCutoff,
                               pAdjustMethod=pAdjustMethod,
                               universe = universe,
                               qvalueCutoff = qvalueCutoff,
                               minGSSize = minGSSize,
                               maxGSSize = maxGSSize,
                               USER_DATA = GO_DATA_drug
      )
      
      if (is.null(res))
        return(res)
      
      res@organism <- get_organism(OrgDb="org.Hs.eg.db")
      res@ontology <- ont
      res@targets <- NULL
      if (ont == "ALL") {
        res <- add_GO_Ontology(res, GO_DATA_drug)
      }
      return(res)
    }
    if(type == "KEGG"){
      species <- organismMapper("hsa")
      KEGG_DATA_drug <- prepare_KEGG_drug(species, "KEGG", keyType="kegg")
      # get all the drugs in the corresponding annotation system as universe
      ext2path <- get("EXTID2PATHID", envir = KEGG_DATA_drug)
      universe = names(ext2path)
      res <- enricher_internal(drugs,
                               pvalueCutoff  = pvalueCutoff,
                               pAdjustMethod = pAdjustMethod,
                               universe      = universe,
                               minGSSize     = minGSSize,
                               maxGSSize     = maxGSSize,
                               qvalueCutoff  = qvalueCutoff,
                               USER_DATA = KEGG_DATA_drug)
      if (is.null(res))
        return(res)
      res@targets <- NULL
      res@ontology <- "KEGG"
      res@organism <- species
      return(res)
    }
}


##' @importFrom AnnotationDbi keys
##' @importFrom AnnotationDbi select
##' @importFrom AnnotationDbi keytypes
##' @importFrom AnnotationDbi toTable
##' @importFrom GO.db GOTERM
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
get_GO_data_drug <- function(OrgDb, ont, keytype) {
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
    
    # # get dtlink in drugbank, lincs and stitch databases
    # ext_path <- system.file("extdata", package="signatureSearch")
    # dtlink_path <- paste0(ext_path,"/dtlink_db_lincs_sti.db")
    # if(file.exists(dtlink_path)){
    #   conn <- dbConnect(SQLite(), dtlink_path)
    # } else {
    #   tryCatch(download.file("http://biocluster.ucr.edu/~yduan004/DOSE2/dtlink_db_lincs_sti.db", dtlink_path, quiet = TRUE), 
    #            error = function(e){
    #              stop("Error happens when downloading signatureSearch/dtlink_db_lincs_sti.db")
    #            }, warning = function(w) {file.remove(dtlink_path)
    #              stop("Error happens when downloading signatureSearch/dtlink_db_lincs_sti.db")})
    #   conn <- dbConnect(SQLite(), dtlink_path)
    # }
    # dtlink <- dbGetQuery(conn, 'SELECT * FROM dtlink')
    # dbDisconnect(conn)
    # 
    # # convert SYMBOL-GOALL in goAnno to drug-GOALL
    # message("left join goAnno and dtlink")
    # goAnno_drug <- left_join(as_tibble(goAnno[,c(2,1,4)]), as_tibble(dtlink), by = c("SYMBOL"="t_gn_sym"))
    # rm(dtlink)
    # goAnno_drug <- as.data.frame(goAnno_drug)[,c("GOALL","ONTOLOGYALL","drug_name")]
    # message("remove na in goAnno_drug")
    # goAnno_drug <- na.omit(goAnno_drug)
    # message("unique GOALL and drug_name in goAnno_drug")
    # goAnno_drug <- goAnno_drug[!duplicated(goAnno_drug[,c("GOALL","drug_name")]),]
    # message("unique is done")
    # # goAnno_drug <- goAnno_drug %>% filter(!is.na(drug_name)) %>% distinct(GOALL, drug_name, .keep_all = TRUE)
    # # goAnno_drug <- as.data.frame(goAnno_drug)
    
    # download goAnno_drug.rds
    ext_path <- system.file("extdata", package="signatureSearch")
    anno_path <- file.path(ext_path, "goAnno_drug.rds")
    if(file.exists(anno_path)){
      goAnno_drug <- readRDS(anno_path)
    } else {
      download.file("http://biocluster.ucr.edu/~yduan004/fea/goAnno_drug.rds", anno_path, quiet = TRUE)
      goAnno_drug <- readRDS(anno_path)
    }
    # "drug_name" in goAnno_drug are all lowercase
    if (ont == "ALL") {
        GO2GENE <- goAnno_drug[,c("GOALL","drug_name")]
    } else {
        GO2GENE <- goAnno_drug[goAnno_drug$ONTOLOGYALL == ont, c("GOALL","drug_name")]
    }

    GO_DATA <- build_Anno(GO2GENE, get_GO2TERM_table())

    goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
    goOnt <- goOnt.df[,2]
    names(goOnt) <- goOnt.df[,1]
    assign("GO2ONT", goOnt, envir=GO_DATA)
    return(GO_DATA)
}

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
prepare_KEGG_drug <- function(species, KEGG_Type="KEGG", keyType="kegg") {
  kegg <- download_KEGG(species, KEGG_Type, keyType)
  
  # get dtlink_entrez
  ext_path <- system.file("extdata", package="signatureSearch")
  dtlink_path <- paste0(ext_path,"/dtlink_db_lincs_sti.db")
  if(file.exists(dtlink_path)){
    conn <- dbConnect(SQLite(), dtlink_path)
  } else {
    tryCatch(download.file("http://biocluster.ucr.edu/~yduan004/DOSE2/dtlink_db_lincs_sti.db", dtlink_path, quiet = TRUE), 
             error = function(e){
               stop("Error happens when downloading DOSE2/dtlink_db_lincs_sti.db")
             }, warning = function(w) {file.remove(dtlink_path)
               stop("Error happens when downloading DOSE2/dtlink_db_lincs_sti.db")})
    conn <- dbConnect(SQLite(), dtlink_path)
  }
  dtlink_entrez <- dbGetQuery(conn, 'SELECT * FROM dtlink_entrez')
  dbDisconnect(conn)
  
  # map entrez id to drug name
  keggpath2entrez <- kegg$KEGGPATHID2EXTID
  keggpath2entrez <- left_join(as_tibble(keggpath2entrez), as_tibble(dtlink_entrez), by = c("to"="ENTREZID")) %>% filter(!is.na(drug_name)) %>% distinct(from, drug_name, .keep_all = TRUE)
  keggpath2entrez <- as.data.frame(keggpath2entrez)[,c("from","drug_name")]
  kegg$KEGGPATHID2EXTID <- keggpath2entrez
  build_Anno(kegg$KEGGPATHID2EXTID,
             kegg$KEGGPATHID2NAME)
}

## get rid of "Undefined global functions or variables" note
cell = drug_name = from = . = NULL
