##' This function can be used to get gene/protein target annotation of query 
##' drugs in specified drug-target annotation resource, such as DrugBank, CLUE 
##' and STITCH. A SQLite database storing drug-target links in the above 
##' three databases can be found at 
##' \code{\link[signatureSearchData]{signatureSearchData}} package. 
##'
##' @title Get gene/protein targets of query drugs
##' @param drugs character vector of drug names
##' @param database drug-target annotation resource. one of 'DrugBank', 'CLUE', 
##' 'STITCH' or 'all'. If 'all', the targets from DrugBank, CLUE 
##' and STITCH databases will be combined.  
##' @return data.frame, one column is query drugs, the other column is
##' target gene symbols.
##' @importFrom RSQLite dbConnect
##' @importFrom RSQLite dbGetQuery
##' @importFrom RSQLite SQLite
##' @importFrom RSQLite dbDisconnect
##' @importFrom signatureSearchData load_sqlite
##' @seealso \code{\link[signatureSearchData]{dtlink_db_clue_sti}}
##' @examples 
##' data(drugs)
##' dt <- get_targets(drugs)
##' @export 

get_targets <- function(drugs, database="all"){
  drugs_orig <- unique(drugs)
  drugs <- unique(tolower(drugs))
  # load dtlink_db_clue_sti.db stored in AnnotationHub
  conn <- load_sqlite("AH69083")
  dtlink_db <- dbGetQuery(conn, 'SELECT * FROM dtlink_db')
  dtlink_clue <- dbGetQuery(conn, 'SELECT * FROM dtlink_clue')
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
      message("No targets found in DrugBank database for ", 
              length(drugs_notar_db),
              " drugs: \n",
              paste(drugs_notar_db, collapse = " / "))
    return(res_db)
  }
  
  # drug targets in LINCS
  dtlist_lincs <- split(dtlink_clue$t_gn_sym, dtlink_clue$drug_name)
  drugs_notar_lincs <- setdiff(drugs, names(dtlist_lincs))
  dtlist_lincs_drugs <- dtlist_lincs[intersect(names(dtlist_lincs), drugs)]
  res_lincs <- list2slash(dtlist_lincs_drugs)
  idx_lincs <- match(res_lincs$drug_name, drugs)
  res_lincs$drug_name = drugs_orig[idx_lincs]
  if(database=="CLUE"){
    if(length(drugs_notar_lincs) > 0)
      message("No targets found in LINCS database for ", 
              length(drugs_notar_lincs), 
              " drugs: \n",
              paste(drugs_notar_lincs, collapse = " / "))
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
      message("No targets found in STITCH database for ", 
              length(drugs_notar_sti), 
              " drugs: \n",
              paste(drugs_notar_sti, collapse = " / "))
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
      message("No targets found in DrugBank/CLUE/STITCH database for ", 
              length(drugs_notar), 
              " drugs: \n",
              paste(drugs_notar, collapse = " / "), "\n")
    return(res)
  }
}

slash2link <- function(slash){
  res <- data.frame()
  for(i in seq_len(nrow(slash))){
    tar <- unlist(strsplit(as.character(slash[i,2]), "; "))
    tmp <- data.frame(drug_name=slash[i,1], t_gn_sym=tar, 
                      stringsAsFactors = FALSE)
    res <- rbind(res, tmp)
  }
  res <- unique(res)
  return(res)
}

list2slash <- function(list){
  if(length(list)==0){
    return(data.frame(drug_name=NULL, t_gn_sym=NULL))
  }
  tar <- vapply(list, function(x) paste0(x, collapse = "; "), 
                FUN.VALUE = character(1))
  drug <- names(list)
  res <- data.frame(drug_name=drug, t_gn_sym=tar, stringsAsFactors = FALSE, 
                    row.names = NULL)
  res <- res[res$t_gn_sym != "", ]
  return(res)
}
