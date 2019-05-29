##' Summarize GESS results on MOA level
##' 
##' Summarize GESS results to get a list of MOA categories ranked by their 
##' connectivity to query signature in cells
##' @param gess_tb table in \code{\link{gessResult}} object
##' @param moa_cats if set as "default", it uses MOA annotations from CLUE 
##' website \url{https://clue.io/}. 
##' Users can customize it as a `list` with names as MOA categories and 
##' contents as drug names.
##' @param cells one of "normal", "cancer" or "all", or a character vector 
##' containing selected cell types. It represents a group of interested 
##' cell types. 
##' \itemize{
##'   \item "all": all the 30 cell types in LINCS database;
##'   \item "normal": normal cell types in LINCS database as one group; 
##'   \item "tumor": tumor cell types in LINCS database as one group;
##' }
##'  
##' @return data.frame 
##' @seealso \code{\link{gessResult}}
##' @examples 
##' db_dir <- system.file("extdata", "sample_db", package = "signatureSearch")
##' sample_db <- loadHDF5SummarizedExperiment(db_dir)
##' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
##' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
##' qsig_fisher <- qSig(query=query_mat, gess_method="Fisher", refdb=sample_db)
##' fisher <- gess_fisher(qSig=qsig_fisher, higher=1, lower=-1)
##' res_moa <- moa_conn(result(fisher))
##' @export
moa_conn <- function(gess_tb, moa_cats="default", cells="normal"){
  if(moa_cats=="default"){
    dir <- system.file("extdata", package = "signatureSearch")
    moa_path <- file.path(dir, "clue_moa_list.rds") 
    moa_list <- readRDS(moa_path)
  }
  if(is(moa_cats, "list")){
    moa_list <- moa_cats
  }
  moa_mrk <- moa_mrk(gess_tb, moa_list, cells=cells)
  return(moa_mrk)
}

##' @importFrom dplyr mutate
moa_mrk <- function(gess_tb, moa_list, cells){
  gess_tb %<>% mutate(rank=seq_len(nrow(gess_tb)))
  if(cells=="all"){
    moa_mrk <- sort(vapply(moa_list, function(x) 
      mean(as.numeric(gess_tb$rank[gess_tb$pert %in% x]), na.rm=TRUE),
      FUN.VALUE = numeric(1)))
    moa_mrk_df <- data.frame(moa=names(moa_mrk), cell="all", 
                             mean_rank=moa_mrk, row.names = NULL, 
                             stringsAsFactors = FALSE)
  }
  else if(cells=="normal"){
    cell_sel <- c("HA1E", "HCC515", "HEK293T", "MCF10A", "NKDBA", "NEU", "NPC", 
                  "FIBRNPC", "ASC", "CD34", "PHH", "SKB")
    gess_tb_sub <- filter(gess_tb, cell %in% cell_sel)
    moa_mrk <- sort(vapply(moa_list, function(x) 
      mean(as.numeric(gess_tb_sub$rank[gess_tb_sub$pert %in% x]), na.rm=TRUE),
      FUN.VALUE = numeric(1)))
    moa_mrk_df <- data.frame(moa=names(moa_mrk), cells="normal", 
                             mean_rank=moa_mrk, row.names = NULL, 
                             stringsAsFactors = FALSE)
  }
  else if(cells=="tumor"){
    cell_sel <- c("A375", "A549", "BT20", "HEPG2", "HL60", "HS578T", "HT29", 
                  "HUH7", "JURKAT", "MCF7", "MDAMB231", "NOMO1", "PC3",
                  "SKBR3", "THP1", "U266", "U937", "VCAP")
    gess_tb_sub <- filter(gess_tb, cell %in% cell_sel)
    moa_mrk <- sort(vapply(moa_list, function(x) 
      mean(as.numeric(gess_tb_sub$rank[gess_tb_sub$pert %in% x]), na.rm=TRUE),
      FUN.VALUE = numeric(1)))
    moa_mrk_df <- data.frame(moa=names(moa_mrk), cells="tumor", 
                             mean_rank=moa_mrk, row.names = NULL, 
                             stringsAsFactors = FALSE)
  } else {
    cell_sel <- cells
    gess_tb_sub <- filter(gess_tb, cell %in% cell_sel)
    moa_mrk <- sort(vapply(moa_list, function(x) 
      mean(as.numeric(gess_tb_sub$rank[gess_tb_sub$pert %in% x]), na.rm=TRUE),
      FUN.VALUE = numeric(1)))
    moa_mrk_df <- data.frame(moa=names(moa_mrk), cells="group", 
                             mean_rank=moa_mrk, row.names = NULL, 
                             stringsAsFactors = FALSE)
  }
  moa_num <- data.frame(moa=names(moa_list), 
                        n_drug=vapply(moa_list, length, FUN.VALUE = integer(1)),
                        stringsAsFactors = FALSE)
  res <- as_tibble(moa_mrk_df) %>% left_join(as_tibble(moa_num), by="moa")
  return(res)
}
