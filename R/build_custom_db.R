#' Build custom reference signature database for GESS 
#' 
#' @title build_custom_db
#' @aliases build_custom_db
#' @param df data.frame, represents genome-wide Gene Expression Profiles (GEPs)
#' (log2FC, z-scores, intensity values, etc.) of compound or genetic treatments.
#' Rownames (Gene SYMBOLs) should be included. The colnames should be of 
#' `(drug)__(cell)__(factor)` format, e.g., `sirolimus__MCF7__trt_cp`.
#' @param dest_dir directory path to store the reference database, which is 
#' the hdf5 backed \code{\link{SummarizedExperiment}} object
#' @return hdf5 backed \code{\link{SummarizedExperiment}} object
#' @seealso \code{\link{saveHDF5SummarizedExperiment}}
#' @importFrom readr read_tsv
#' @importFrom stats rnorm
#' @examples 
#' # Generate a data.frame 
#' df <- data.frame(sirolimus__MCF7__trt_cp=rnorm(1000),
#'                  vorinostat__SKB__trt_cp=rnorm(1000))
#' data(targetList)
#' rownames(df) = names(targetList)
#' dest_dir = tempfile("h5_df_")
#' # The following step takes some time, so do not run as example
#' # se = build_custom_db(df, dest_dir)
#' list.files(dest_dir)
#' unlink(dest_dir, recursive=TRUE)
#' @export

build_custom_db <- function(df, dest_dir){
  se <- SummarizedExperiment(assays = list(score=as.matrix(df)))
  se <- saveHDF5SummarizedExperiment(se, dest_dir, replace=TRUE)
  return(se)
}
