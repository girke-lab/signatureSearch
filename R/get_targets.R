##' This function returns for a set of query drug names/ids the corresponding
##' target gene/protein ids. The required drug-target annotations are from
##' DrugBank, CLUE and STITCH. An SQLite database storing these drug-target
##' interactions based on the above three annotation resources is available in 
##' the \code{\link[signatureSearchData]{signatureSearchData}} package. 
##'
##' @title Target Gene/Protein IDs for Query Drugs
##' @param drugs character vector of drug names
##' @param database drug-target annotation resource; one of 'DrugBank', 'CLUE',
##' 'STITCH' or 'all'. If 'all', the targets from DrugBank, CLUE 
##' and STITCH databases will be combined.  
##' @param verbose TRUE or FALSE, whether to print messages
##' @param collapse TRUE or FALSE. If setting as FALSE, the result is a 
##' data.frame of target gene symbols of each drug. If TRUE, the targets are 
##' collapsed to be a character vector with duplications if different drugs 
##' have the same targets.
##' @return data.frame, one column contains the query drug names and the other 
##' target gene symbols if \code{collapse} is FALSE. A character vetor if
##' \code{collapse} is TRUE.
##' @importFrom RSQLite dbConnect
##' @importFrom RSQLite dbGetQuery
##' @importFrom RSQLite SQLite
##' @importFrom RSQLite dbDisconnect
##' @seealso \code{\link[signatureSearchData]{dtlink_db_clue_sti}}
##' @examples 
##' data(drugs10)
##' dt <- get_targets(drugs10)
##' @export 

get_targets <- function(drugs, database="all", verbose=TRUE, collapse=FALSE){
  drugs_orig <- unique(drugs)
  drugs <- unique(tolower(drugs))
  # load dtlink_db_clue_sti.db stored in AnnotationHub
  conn <- load_sqlite("EH3228")
  dtlink_db <- dbGetQuery(conn, 'SELECT * FROM dtlink_db')
  dtlink_clue <- dbGetQuery(conn, 'SELECT * FROM dtlink_clue')
  dtlink_sti <- dbGetQuery(conn, 'SELECT * FROM dtlink_sti')
  dtlink <- dbGetQuery(conn, 'SELECT * FROM dtlink')
  dbDisconnect(conn)
  
  # drug targets in DrugBank
  if(database=="DrugBank"){
      dtlist_db <- split(dtlink_db$t_gn_sym, dtlink_db$drug_name)
      drugs_notar_db <- setdiff(drugs, names(dtlist_db))
      dtlist_db_drugs <- dtlist_db[intersect(drugs, names(dtlist_db))]
      res <- list2slash(dtlist_db_drugs)
      idx_db <- match(res$drug_name, drugs)
      res$drug_name <- drugs_orig[idx_db]
      if(length(drugs_notar_db) > 0 & verbose){
          message("No targets found in DrugBank database for ", 
                  length(drugs_notar_db),
                  " drugs: \n",
                  paste(drugs_notar_db, collapse = " / "))
      }
  }
  # drug targets in CLUE
  if(database=="CLUE"){
      dtlist_clue <- split(dtlink_clue$t_gn_sym, dtlink_clue$drug_name)
      drugs_notar_clue <- setdiff(drugs, names(dtlist_clue))
      dtlist_clue_drugs <- dtlist_clue[intersect(drugs, names(dtlist_clue))]
      res <- list2slash(dtlist_clue_drugs)
      idx_clue <- match(res$drug_name, drugs)
      res$drug_name <- drugs_orig[idx_clue]
      
      if(length(drugs_notar_clue) > 0 & verbose){
          message("No targets found in CLUE database for ", 
                  length(drugs_notar_clue), 
                  " drugs: \n",
                  paste(drugs_notar_clue, collapse = " / "))
      }
  }
  
  # drug targets in STITCH
  if(database=="STITCH"){
      dtlist_sti <- split(dtlink_sti$t_gn_sym, dtlink_sti$drug_name)
      drugs_notar_sti <- setdiff(drugs, names(dtlist_sti))
      dtlist_sti_drugs <- dtlist_sti[intersect(drugs, names(dtlist_sti))]
      res <- list2slash(dtlist_sti_drugs)
      idx_sti <- match(res$drug_name, drugs)
      res$drug_name <- drugs_orig[idx_sti]
    
      if(length(drugs_notar_sti) > 0 & verbose){
          message("No targets found in STITCH database for ", 
                  length(drugs_notar_sti), 
                  " drugs: \n",
                  paste(drugs_notar_sti, collapse = " / "))
      }
  }
  
  # drug targets in all database
  if(database=="all"){
      dtlist <- split(dtlink$t_gn_sym, dtlink$drug_name)
      drugs_notar <- setdiff(drugs, names(dtlist))
      dtlist_drugs <- dtlist[intersect(drugs, names(dtlist))]
      res <- list2slash(dtlist_drugs)
      idx <- match(res$drug_name, drugs)
      res$drug_name <- drugs_orig[idx]
  
      if(length(drugs_notar) > 0 & verbose){
          message("No targets found in DrugBank/CLUE/STITCH database for ", 
                  length(drugs_notar), 
                  " drugs: \n",
                  paste(drugs_notar, collapse = " / "), "\n")
      }
  }
  
  if(collapse){
      gnset <- na.omit(unlist(lapply(res$t_gn_sym, function(i) 
          unlist(strsplit(as.character(i), split = "; ")))))
      return(gnset)
  } else {
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
