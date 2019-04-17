#' Build custom reference signature database for GESS methods
#' 
#' It write a data frame containing GEPs of treatments into an HDF5 file as
#' reference database
#' 
#' @title build_custom_db
#' @aliases build_custom_db
#' @param df data.frame, represents genome-wide Gene Expression Profiles (GEPs)
#' of compound or genetic treatments in cells. The GEPs can be log2 fold change,
#' z-scores \emph{etc.} from differential expression analysis, or gene 
#' expression intensity values from Affymetrix Chips, or read counts from 
#' RNA-Seq data. 
#' 
#' Rownames (gene Entrez id) should be included. The colnames should be of 
#' `(drug)__(cell)__(factor)` format, e.g., `sirolimus__MCF7__trt_cp`.
#' It can be generalized to any type of treatment, for example, gene knockdown
#' or over expression by setting `drug` as gene name, `factor` as 'gene_ko' or
#' other words, `cell` as treatment cell types. So, one example for 
#' generalization format could be `P53__MCF7__gene_ko`.  
#' @param h5file character(1), path to the destination hdf5 file
#' @return HDF5 file
#' @importFrom readr read_tsv
#' @importFrom stats rnorm
#' @importFrom signatureSearchData matrix2h5
#' @examples 
#' # Generate a data.frame 
#' df <- data.frame(sirolimus__MCF7__trt_cp=rnorm(1000),
#'                  vorinostat__SKB__trt_cp=rnorm(1000))
#' data(targetList)
#' rownames(df) = names(targetList)
#' h5file = tempfile(fileext=".h5")
#' build_custom_db(df, h5file)
#' library(signatureSearchData)
#' tmp <- readHDF5chunk(h5file, colindex=1:2)
#' @export

build_custom_db <- function(df, h5file){
    mat <- as.matrix(df)
    matrix2h5(mat, h5file, overwrite=TRUE)
}
