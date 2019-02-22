##' get targets in DrugBank/LINCS/STITCH for a given list of drugs
##'
##' @title get drug targets in DrugBank/LINCS/STITCH
##' @param drugs character vector, a list of drug names
##' @param database one of "DrugBank", "LINCS", "STITCH" or "all"  
##' @return data.frame of drugs and target symbols
##' @importFrom drugbankR queryDB
##' @importFrom RSQLite dbConnect
##' @importFrom RSQLite dbGetQuery
##' @importFrom utils download.file
##' @importFrom RSQLite SQLite
##' @importFrom RSQLite dbDisconnect
##' @export 

get_targets <- function(drugs, database="all"){
  drugs_orig <- drugs
  drugs <- unique(tolower(drugs))
  ext_path <- system.file("extdata", package="signatureSearch")
  dtlink_path <- paste0(ext_path,"/dtlink_db_lincs_sti.db")
  if(file.exists(dtlink_path)){
    conn <- dbConnect(SQLite(), dtlink_path)
  } else {
    tryCatch(download.file("http://biocluster.ucr.edu/~yduan004/DOSE2/dtlink_db_lincs_sti.db", dtlink_path, quiet = TRUE), 
             error = function(e){
               stop("Error happens when downloading signatureSearch/dtlink_db_lincs_sti.db")
             }, warning = function(w) {file.remove(dtlink_path)
               stop("Error happens when downloading signatureSearch/dtlink_db_lincs_sti.db")})
    conn <- dbConnect(SQLite(), dtlink_path)
  }
  dtlink_db <- dbGetQuery(conn, 'SELECT * FROM dtlink_db')
  dtlink_lincs <- dbGetQuery(conn, 'SELECT * FROM dtlink_lincs')
  dtlink_sti <- dbGetQuery(conn, 'SELECT * FROM dtlink_sti')
  dtlink <- dbGetQuery(conn, 'SELECT * FROM dtlink')
  dbDisconnect(conn)
  
  # drug targets in DrugBank
  dtlist_db <- split(dtlink_db$t_gn_sym, dtlink_db$drug_name)
  drugs_notar_db <- setdiff(drugs, names(dtlist_db))
  dtlist_db_drugs <- dtlist_db[intersect(names(dtlist_db), drugs)]
  res_db <- list2slash(dtlist_db_drugs)
  idx_db <- match(res_db$drug_name, drugs)
  res_db$drug_name = drugs_orig[idx_db]
  if(database=="DrugBank"){
    if(length(drugs_notar_db) > 0)
      message("No targets found in DrugBank database for ", length(drugs_notar_db)," drugs: \n",paste(drugs_notar_db, collapse = " / "))
    return(res_db)
  }
  
  # drug targets in LINCS
  dtlist_lincs <- split(dtlink_lincs$t_gn_sym, dtlink_lincs$drug_name)
  drugs_notar_lincs <- setdiff(drugs, names(dtlist_lincs))
  dtlist_lincs_drugs <- dtlist_lincs[intersect(names(dtlist_lincs), drugs)]
  res_lincs <- list2slash(dtlist_lincs_drugs)
  idx_lincs <- match(res_lincs$drug_name, drugs)
  res_lincs$drug_name = drugs_orig[idx_lincs]
  if(database=="LINCS"){
    if(length(drugs_notar_lincs) > 0)
      message("No targets found in LINCS database for ", length(drugs_notar_lincs), " drugs: \n",paste(drugs_notar_lincs, collapse = " / "))
    return(res_lincs)
  }
  
  # drug targets in STITCH
  dtlist_sti <- split(dtlink_sti$t_gn_sym, dtlink_sti$drug_name)
  drugs_notar_sti <- setdiff(drugs, names(dtlist_sti))
  dtlist_sti_drugs <- dtlist_sti[intersect(names(dtlist_sti), drugs)]
  res_sti <- list2slash(dtlist_sti_drugs)
  idx_sti <- match(res_sti$drug_name, drugs)
  res_sti$drug_name = drugs_orig[idx_sti]
  if(database=="STITCH"){
    if(length(drugs_notar_sti) > 0)
      message("No targets found in STITCH database for ", length(drugs_notar_sti), " drugs: \n",paste(drugs_notar_sti, collapse = " / "))
    return(res_sti)
  }
  
  # drug targets in all database
  dtlist <- split(dtlink$t_gn_sym, dtlink$drug_name)
  drugs_notar <- setdiff(drugs, names(dtlist))
  dtlist_drugs <- dtlist[intersect(names(dtlist), drugs)]
  res <- list2slash(dtlist_drugs)
  idx <- match(res$drug_name, drugs)
  res$drug_name = drugs_orig[idx]
  if(database=="all"){
    if(length(drugs_notar) > 0)
      message("No targets found in DrugBank/LINCS/STITCH database for ", length(drugs_notar), " drugs: \n",paste(drugs_notar, collapse = " / "), "\n")
    return(res)
  }
}

slash2link <- function(slash){
  res <- data.frame()
  for(i in 1:dim(slash)[1]){
    tar <- unlist(strsplit(as.character(slash[i,2]), "; "))
    tmp <- data.frame(drug_name=slash[i,1], t_gn_sym=tar, stringsAsFactors = FALSE)
    res <- rbind(res, tmp)
  }
  res <- unique(res)
  return(res)
}

list2slash <- function(list){
  tar <- sapply(list, function(x) paste0(x, collapse = "; "))
  drug <- names(list)
  res <- data.frame(drug_name=drug, t_gn_sym=tar, stringsAsFactors = FALSE, row.names = NULL)
  res <- res[res$t_gn_sym != "",]
  return(res)
}