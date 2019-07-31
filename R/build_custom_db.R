#' Build custom reference signature database for GESS methods
#' 
#' It stores a data.frame or matrix containing genome-wide gene expression 
#' data (e.g. for drug, disease or genetic perturbations) into an HDF5 file as
#' reference database. The gene expression data can be most types of the 
#' pre-processed gene expression values, such as gene expression intensity 
#' values (or counts for RNA-Seq), log2 fold changes (LFC), z-scores or 
#' p-values obtained from DE analysis.
#' 
#' @title build_custom_db
#' @aliases build_custom_db
#' @param df data.frame or matrix, represents genome-wide GESs
#' of compound or genetic treatments in cells.
#' 
#' Rownames representing gene IDs (e.g. Entrez ids) should be included. 
#' The colnames should be of `(drug)__(cell)__(factor)` format, e.g., 
#' `sirolimus__MCF7__trt_cp`. It can be generalized to any type of treatment, 
#' for example, gene knockdown or over expression by setting `drug` as gene 
#' name, `factor` as 'ko' or other words, `cell` as treatment cell types. So, 
#' one example for generalization format could be `P53__MCF7__ko`.  
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
    # Validity check of df
    if(is.null(rownames(mat))){
        stop(paste("The input data.frame or matrix should have rownames",
        "slot as gene IDs!"))
    }
    if(is.null(colnames(mat))){
        stop(paste("The input data.frame or matrix should have colnames",
                   "slot in `(drug)__(cell)__(factor)` format!"))
    }
    if(length(strsplit(colnames(mat)[1], split="__")[[1]]) != 3){
        stop(paste("The input data.frame or matrix should have colnames",
                   "slot in `(drug)__(cell)__(factor)` format!"))
    }
    matrix2h5(mat, h5file, overwrite=TRUE)
}
