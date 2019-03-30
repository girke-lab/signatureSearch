##' hyperG method for DSEA
##' 
##' The hypergeometric test is used to do enrichment analysis on a drug set 
##' after mapping drugs to functional categories via drug-target links in
##' DrugBank, CLUE and STITCH databases
##' 
##' The drug sets can also be directly used for functional enrichment testing 
##' by changing the mappings in the reference database from target-to-functional
##' category mappings to drug-to-functional category mappings. 
##' The latter can be generated based on the drug-target information provided 
##' by DrugBank or related databases. As a result, one can perform the FEA on 
##' ranked drug lists directly. This DSEA approach has the advantage that the 
##' drugs in the query test sets are usually unique allowing to use them 
##' without any changes for functional enrichment methods. 
##' 
##' So the hypergeometric test is directly used when the query is a vector of 
##' drugs, e.g., top ranking drugs from GESS result since the functional
##' categories also contains drug sets after gene to drug mappings.
##'
##' @param drugs query drug set used to do DSEA. 
##' Can be top ranking drugs in GESS result. 
##' @param type one of "GO" or "KEGG"
##' @param ont One of "MF", "BP", and "CC" or "ALL".
##' @param pvalueCutoff Cutoff value of p value.
##' @param pAdjustMethod p value adjustment method,
##' one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param qvalueCutoff q value cutoff
##' @param minGSSize minimal size of drug sets annotated by ontology term 
##' after drug to functional category mappings.
##' @param maxGSSize maximal size of drug sets annotated for testing
##' @return \code{\link{feaResult}} object, 
##' represents enriched functional categories.
##' @import org.Hs.eg.db
##' @importFrom magrittr %<>%
##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importMethodsFrom AnnotationDbi mget
##' @seealso \code{\link{feaResult}}, \code{\link{fea}},
##'          \code{\link[signatureSearch_data]{dtlink_db_clue_sti.db}}
##' @examples 
##' data(drugs)
##' # GO annotation system
##' hyperG_res <- dsea_hyperG(drugs = drugs, type = "GO", ont="MF")
##' result(hyperG_res)
##' ## KEGG annotation system
##' hyperG_k_res <- dsea_hyperG(drugs = drugs, type = "KEGG", 
##'                             pvalueCutoff = 1, qvalueCutoff = 1, 
##'                             minGSSize = 10, maxGSSize = 2000)
##' result(hyperG_k_res) 
##' @export
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
    
    # GO_DATA_drug <- get_GO_data_drug(OrgDb = "org.Hs.eg.db", 
    #                                  ont, keytype="SYMBOL")
    # download GO_DATA_drug.rds and save it to cache to save time
    fl <- download_data_file(url=
    "http://biocluster.ucr.edu/~yduan004/signatureSearch_data/GO_DATA_drug.rds",
                             rname="GO_DATA_drug")
    GO_DATA_drug <- readRDS(fl)
    
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
                             USER_DATA = GO_DATA_drug)
    
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
  
  # download goAnno_drug.rds and save it to cache
  fl <- download_data_file(url=
    "http://biocluster.ucr.edu/~yduan004/signatureSearch_data/goAnno_drug.rds",
                           rname="goAnno_drug")
  goAnno_drug <- readRDS(fl) # "drug_name" in goAnno_drug are all lowercase
  
  if (ont == "ALL") {
    GO2GENE <- goAnno_drug[,c("GOALL","drug_name")]
  } else {
    GO2GENE <- goAnno_drug[goAnno_drug$ONTOLOGYALL == ont, 
                           c("GOALL","drug_name")]
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
##' @importFrom clusterProfiler download_KEGG
prepare_KEGG_drug <- function(species, KEGG_Type="KEGG", keyType="kegg") {
  kegg <- clusterProfiler::download_KEGG(species, KEGG_Type, keyType)
  # get dtlink_entrez
  fl <- download_data_file(url=paste0("http://biocluster.ucr.edu/~yduan004/",
                                  "signatureSearch_data/dtlink_db_clue_sti.db"),
                           rname="dtlink")
  conn <- dbConnect(SQLite(), fl)
  dtlink_entrez <- dbGetQuery(conn, 'SELECT * FROM dtlink_entrez')
  dbDisconnect(conn)
  
  # map entrez id to drug name
  keggpath2entrez <- kegg$KEGGPATHID2EXTID
  keggpath2entrez <- left_join(as_tibble(keggpath2entrez), 
                               as_tibble(dtlink_entrez), 
                               by = c("to"="ENTREZID")) %>% 
    filter(!is.na(drug_name)) %>% distinct(from, drug_name, .keep_all = TRUE)
  keggpath2entrez <- as.data.frame(keggpath2entrez)[,c("from","drug_name")]
  kegg$KEGGPATHID2EXTID <- keggpath2entrez
  build_Anno(kegg$KEGGPATHID2EXTID, kegg$KEGGPATHID2NAME)
}
## get rid of "Undefined global functions or variables" note
cell = drug_name = from = . = NULL
