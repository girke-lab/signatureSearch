#' Build custom reference signature database for GESS 
#' 
#' @title build_custom_db
#' @aliases build_custom_db
#' @param df data.frame, represents genome-wide Gene Expression Profiles (GEPs)
#' of compound or genetic treatments in cells. The GEPs can be log2 fold change,
#' z-scores \emph{etc.} from differential expression analysis, or gene 
#' expression intensity values from Affymetrix Chips, or read counts from 
#' RNA-Seq data. 
#' 
#' Rownames (Gene SYMBOLs) should be included. The colnames should be of 
#' `(drug)__(cell)__(factor)` format, e.g., `sirolimus__MCF7__trt_cp`.
#' It can be generalized to any type of treatment, for example, gene knockdown
#' or overexpression by setting `drug` as gene symbol, `factor` as 'gene_ko' or 
#' other words, `cell` as treatment cell types. So, one example for 
#' generalization format could be `P53__MCF7__gene_ko`.  
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
