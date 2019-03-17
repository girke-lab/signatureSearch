#' Build custom reference signature database for GESS 
#' 
#' @title build_custom_db
#' @param df data.frame, represents genoime-wide GEPs (log2FC, z-scores, 
#' intensity values, etc.) of compound or genetic treatments. 
#' Rownames (Gene IDs) should be included. The colnames should be of 
#' `(drug)__(cell)__(factor)` format, e.g., `sirolimus__MCF7__trt_cp`.
#' @param dest_path directory path to store the reference database, which is 
#' the hdf5 backed \code{SummarizedExperiment} object
#' @return hdf5 backed \code{SummarizedExperiment} object
#' @importFrom readr read_tsv
#' @export

build_custom_db <- function(df, dest_path){
  df <- read.delim(file = file, sep = "\t", row.names = 1, check.names = FALSE)
  se <- SummarizedExperiment(assays = list(score=as.matrix(df)))
  se <- saveHDF5SummarizedExperiment(se, dest_path)
  return(se)
}
