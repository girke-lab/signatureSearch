##' This function returns for a set of query drug names/ids the corresponding
##' target gene/protein ids. The required drug-target annotations are from
##' DrugBank, CLUE and STITCH. An SQLite database storing these drug-target
##' interactions based on the above three annotation resources is available in 
##' the \code{\link[signatureSearchData]{signatureSearchData}} package. 
##'
##' @title Target Gene/Protein IDs for Query Drugs
##' @param drugs character vector of drug names
##' @param database drug-target annotation resource; A character vector of any
##' combination of 'DrugBank', 'CLUE', STITCH' or 'all'. The target set from
##' the selected resources will be combined. If 'all' is contained in the 
##' character vector, target sets from all of the annotation databases
##' (DrugBank, CLUE and STITCH) will be combined.
##' @param verbose TRUE or FALSE, whether to print messages
##' @param output one of "df", "list" or "vector". If setting as "df", the result 
##' is in a data.frame format containing target gene symbols separated by semicolon
##' for each drug. If setting as "list", the result is a list of targets for each
##' query drug. If setting as "vector", the result is a character vector of the 
##' target set that are collapsed with duplications if different drugs 
##' have the same targets.
##' @return drug-target annotation in a format defined by the \code{output} argument.
##' @importFrom RSQLite dbConnect
##' @importFrom RSQLite dbGetQuery
##' @importFrom RSQLite SQLite
##' @importFrom RSQLite dbDisconnect
##' @seealso \code{\link[signatureSearchData]{dtlink_db_clue_sti}}
##' @examples 
##' data(drugs10)
##' dt <- get_targets(drugs10)
##' @export 

get_targets <- function(drugs, database="all", verbose=TRUE, output="df"){
  drugs_orig <- unique(drugs)
  drugs <- unique(tolower(drugs))
  # load dtlink_db_clue_sti.db stored in AnnotationHub
  conn <- load_sqlite("EH3228")
  dtlink_db <- dbGetQuery(conn, 'SELECT * FROM dtlink_db')
  dtlink_clue <- dbGetQuery(conn, 'SELECT * FROM dtlink_clue')
  dtlink_sti <- dbGetQuery(conn, 'SELECT * FROM dtlink_sti')
  dtlink <- dbGetQuery(conn, 'SELECT * FROM dtlink')
  dbDisconnect(conn)
  
  dtl <- lapply(database, function(x){
      if(tolower(x)=="drugbank"){
          return(dtlink_db)
      }
      if(tolower(x)=="stitch"){
          return(dtlink_sti)
      }
      if(tolower(x)=="clue"){
          return(dtlink_clue)
      }
      if(tolower(x)=="all"){
          return(dtlink)
      }
  })
  dt <- unique(do.call(rbind, dtl))
  
  dtlist <- split(dt$t_gn_sym, dt$drug_name)
  drugs_notar <- setdiff(drugs, names(dtlist))
  if(length(drugs_notar) > 0 & verbose){
      message("No targets found in ", paste(database, collapse="/"), 
              " databases for ", length(drugs_notar), " drugs: \n",
              paste(drugs_notar, collapse = " / "), "\n")
  }
  dtlist_drugs <- dtlist[intersect(drugs, names(dtlist))]
  names(dtlist_drugs) <- drugs_orig[match(names(dtlist_drugs), drugs)]
  if(output=="list"){
      return(dtlist_drugs)
  }
  if(output=="vector"){
      return(unlist(dtlist_drugs, use.names=FALSE))
  }
  if(output=="df"){
      res <- list2slash(dtlist_drugs)
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
