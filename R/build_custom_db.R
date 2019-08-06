#' Build custom reference signature database for GESS methods
#' 
#' The perturbation-based gene expression data, here provided as data.frame or 
#' matrix, will be stored in an HDF5 file. The latter can be used as reference 
#' database by compatible GESS methods of signatureSearch. Various types of 
#' pre-processed gene expression data can be used here, such as normalized 
#' gene expression intensities (or counts for RNA-Seq); log2 fold changes (LFC), 
#' Z-scores or p-values obtained from analysis routines of differentially expressed
#' genes (DEGs).
#' 
#' @title build_custom_db
#' @aliases build_custom_db
#' @param df data.frame or matrix containing genome-wide or close to
#' genome-wide GESs of perturbation experiments.
#' 
#' The row name slots are expected to contain gene or transcript IDs 
#' (e.g. Entrez ids), while the column names are expected to have this structure:
#' `(drug)__(cell)__(factor)`, e.g. `sirolimus__MCF7__trt_cp`. This format is 
#' flexible enough to encode most perturbation types of biological samples. For 
#' example, gene knockdown or over expression treatments can be specified by 
#' assigning the ID of the affected gene to `drug`, and `ko` or `ov` to `factor`, 
#' respectively. An example for a knockdown treatment would look like this: 
#' `P53__MCF7__ko`.
#' 
#' @param h5file character vector of length 1 containing the path to the destination 
#'         hdf5 file
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
