sep_pcf <- function(res){
    new <- as.data.frame(t(vapply(seq_len(nrow(res)), function(i)
        unlist(strsplit(as.character(res$set[i]), "__")), 
        FUN.VALUE=character(3))), stringsAsFactors=FALSE)
    colnames(new) = c("pert", "cell", "type")
    res <- cbind(new, res[,-1])
    return(res)
}

#' @import rhdf5
readHDF5mat <- function(h5file, colindex=seq_len(10)) {
    m <- h5read(h5file, "assay", index=list(NULL, colindex))
    mycol <- h5read(h5file, "colnames", index=list(colindex, 1))
    myrow <- h5read(h5file, "rownames")
    h5closeAll()
    rownames(m) <- as.character(myrow[,1])
    colnames(m) <- as.character(mycol[,1])
    return(m)
}

getH5dim <- function(h5file){
    mat_dim <- h5ls(h5file)$dim[1]
    mat_ncol <- as.numeric(gsub(".* x ","", mat_dim))
    mat_nrow <- as.numeric(gsub(" x .*","", mat_dim))
    return(c(mat_nrow, mat_ncol))
}

#' Write Matrix to HDF5 file
#' 
#' Function writes matrix object to an HDF5 file.
#' @param matrix matrix to be written to HDF5 file, row and column name slots 
#' need to be populated
#' @param h5file character(1), path to the hdf5 destination file
#' @param overwrite TRUE or FALSE, whether to overwrite or append
#' matrix to an existing 'h5file'
#' @return HDF5 file containing exported matrix
#' @examples 
#' mat <- matrix(rnorm(12), nrow=3, dimnames=list(
#'               paste0("r",1:3), paste0("c",1:4)))
#' h5file <- tempfile(fileext=".h5")
#' matrix2h5(matrix=mat, h5file=h5file, overwrite=TRUE)
#' @export
matrix2h5 <- function(matrix, h5file, overwrite=TRUE){
    if(file.exists(h5file)){
        if(isTRUE(overwrite)){
            createEmptyH5(h5file, delete_existing=TRUE)
        }
    } else {
        createEmptyH5(h5file, delete_existing=FALSE)
    }
    append2H5(matrix, h5file)
}

#' Create Empty HDF5 File 
#' 
#' This function can be used to create an empty HDF5 file where the user defines
#' the file path and compression level. The empty HDF5 file has under its root
#' group three data slots named 'assay', 'colnames' and 'rownames' for storing a
#' \code{numeric matrix} along with its column names (\code{character}) and row 
#' names (\code{character}), respectively.
#' 
#' @param h5file character(1), path to the HDF5 file to be created
#' @param delete_existing logical, whether to delete an existing HDF5 file with 
#' identical path
#' @param level The compression level used, here given as integer value between 
#' 0 (no compression) and 9 (highest and slowest compression).
#' @return empty HDF5 file
#' @examples
#' tmp_file <- tempfile(fileext=".h5")
#' createEmptyH5(tmp_file, level=6)
#' @export
createEmptyH5 <- function(h5file, delete_existing=FALSE, level=6) {
    if(delete_existing==TRUE) unlink(h5file)
    h5createFile(file=h5file)
    h5createDataset(h5file, "assay", c(0,0), c(H5Sunlimited(), H5Sunlimited()), 
                    chunk=c(12328,1), level=level)
    h5createDataset(h5file, "colnames", c(0,1), c(H5Sunlimited(), 1), 
                    storage.mode='character', size=1000, level=level)
    h5createDataset(h5file, "rownames", c(0,1), c(H5Sunlimited(), 1), 
                    storage.mode='character', size=100, level=level)
}

#' Append Matrix to HDF5 File
#' 
#' Function to write matrix data to an existing HDF5 file. If the file contains 
#' already matrix data then both need to have the same number of rows. The 
#' append will be column-wise.
#' @param x matrix object to write to an HDF5 file. If the HDF5 file is not 
#' empty, the exported matrix data needs to have the same number rows as the 
#' matrix stored in the HDF5 file, and will be appended column-wise to the 
#' existing one.
#' @param h5file character(1), path to existing HDF5 file that can be empty or
#' contain matrix data
#' @param printstatus logical, whether to print status
#' @return HDF5 file storing exported matrix
#' @examples 
#' mat <- matrix(1:12, nrow=3)
#' rownames(mat) <- paste0("r", 1:3); colnames(mat) <- paste0("c", 1:4)
#' tmp_file <- tempfile(fileext=".h5")
#' createEmptyH5(tmp_file)
#' append2H5(mat, tmp_file)
#' rhdf5::h5ls(tmp_file)
#' @export
append2H5 <- function(x, h5file, printstatus=TRUE) {
    status <- h5ls(h5file)[c("name", "dim")]
    rowstatus <- as.numeric(gsub(" x \\d{1,}$", "", 
                                 status[status$name=="assay", "dim"]))
    colstatus <- as.numeric(gsub("^\\d{1,} x ", "", 
                                 status[status$name=="assay", "dim"]))
    nrows <- nrow(x) 
    ncols <- colstatus + ncol(x)
    h5set_extent(h5file, "assay", c(nrows, ncols))
    h5write(x, h5file, "assay", index=list(seq_len(nrows), (colstatus+1):ncols))
    h5set_extent(h5file, "colnames", c(ncols,1))
    h5write(colnames(x), h5file, "colnames", index=list((colstatus+1):ncols, 1))
    if(any(duplicated(h5read(h5file, "colnames")[,1]))) 
        warning("Column names contain duplicates!")
    h5set_extent(h5file, "rownames", c(nrows,1))
    h5write(rownames(x), h5file, "rownames", index=list(seq_len(nrows), 1))
    if(any(duplicated(h5read(h5file, "rownames")[,1]))) 
        warning("Row names contain duplicates!")
    if(printstatus==TRUE) h5ls(h5file, all=TRUE)[c("dim", "maxdim")]
    h5closeAll()
}

#' Read matrix-like data from large gctx file in chunks and write result back 
#' to an HDF5 file.
#' @title Convert GCTX to HDF5 File
#' @param gctx character(1), path to gctx file from LINCS
#' @param cid character or integer vector referencing the
#' columns of the matrix to include
#' @param new_cid character vector of the same length as cid, assigning new
#' column names to matrix
#' @param h5file character(1), path of the hdf5 destination file
#' @param chunksize number of columns to import in each iteration to limit 
#' memory usage
#' @param overwrite TRUE or FALSE, whether to overwrite or to append to 
#' existing 'h5file'
#' @return HDF5 file
#' @import rhdf5
#' @examples 
#' gctx <- system.file("extdata", "test_sample_n2x12328.gctx", 
#'         package="signatureSearch")
#' h5file <- tempfile(fileext=".h5")
#' gctx2h5(gctx, cid=1:2, 
#'         new_cid=c('sirolimus__MCF7__trt_cp', 'vorinostat__SKB__trt_cp'), 
#'         h5file=h5file, chunksize=5000, overwrite=TRUE)
#' @export
#' 
gctx2h5 <- function(gctx, cid, new_cid=cid, h5file, chunksize=5000, 
                    overwrite=TRUE){
    cid_list <- suppressWarnings(
        split(cid, rep(seq_len(ceiling(length(cid)/chunksize)), 
                       each=chunksize)))
    new_cid_list <- suppressWarnings(
        split(new_cid, rep(seq_len(ceiling(length(new_cid)/chunksize)), 
                           each=chunksize)))
    if(file.exists(h5file)){
        if(isTRUE(overwrite)){
            createEmptyH5(h5file, delete_existing=TRUE)
        }
    } else {
        createEmptyH5(h5file, delete_existing=FALSE)
    }
    lapply(seq_along(cid_list), function(i){
        mat <- parse_gctx(gctx, cid=cid_list[[i]], matrix_only=TRUE)
        mat <- mat@mat
        colnames(mat) <- new_cid_list[[i]]
        append2H5(x=mat, h5file, printstatus=FALSE)
    })    
    h5ls(h5file)
}

#' Import HDF5 Data into SummarizedExperiment Object
#' 
#' Imports user-definable subsets of matrix data from an HDF5 file into a 
#' \code{SummarizedExperiment} object. The
#' corresponding HDF5 file is expected to have three data components named
#' 'assay', 'colnames' and 'rownames' containing the numeric values, column
#' names and row names of a matrix, respectively.
#' 
#' @param h5file character(1), path to HDF5 file
#' @param colindex integer vector, position index of the matrix columns to be 
#' imported
#' @param colnames character vector, names of the columns of the matrix to be 
#' imported. If 'colnames' is set, 'colindex' will be ignored.
#' @return \code{SummarizedExperiment} object
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @seealso 
#' \code{SummarizedExperiment}
#' @examples
#' gctx <- system.file("extdata", "test_sample_n2x12328.gctx", 
#'                     package="signatureSearch")
#' h5file <- tempfile(fileext=".h5")
#' gctx2h5(gctx, cid=1:2, 
#'         new_cid=c('sirolimus__MCF7__trt_cp', 'vorinostat__SKB__trt_cp'), 
#'         h5file=h5file, chunksize=5000, overwrite=TRUE)
#' se <- readHDF5chunk(h5file, colindex=1:2)
#' @export
#' 
readHDF5chunk <- function(h5file, colindex=seq_len(10), colnames=NULL) {
    if(! is.null(colnames)){
        all_trts <- h5read(h5file, "colnames", drop=TRUE)
        colindex2 <- which(all_trts %in% colnames)
        m <- h5read(h5file, "assay", index=list(NULL, colindex2))
        colindex <- colindex2
    } else {
        m <- h5read(h5file, "assay", index=list(NULL, colindex))
    }
    mycol <- h5read(h5file, "colnames", index=list(colindex, 1))
    myrow <- h5read(h5file, "rownames")
    rownames(m) <- as.character(myrow[,1])
    colnames(m) <- as.character(mycol[,1])
    if(! is.null(colnames)){
        m = m[,colnames, drop=FALSE]
    }
    se <- SummarizedExperiment(assays=list(score=m))
    h5closeAll()
    return(se)
}

#' @importFrom AnnotationHub AnnotationHub
ah <- NULL
.onLoad <- function(...) {
    ah <<- suppressMessages(AnnotationHub())
}

determine_refdb <- function(refdb){
    if(refdb=="cmap"){
        return(ah[["AH69090"]])
    }
    if(refdb=="cmap_expr"){
        return(ah[["AH69091"]])
    }
    if(refdb=="lincs"){
        return(ah[["AH69092"]])
    }
    if(refdb=="lincs_expr"){
        return(ah[["AH69093"]])
    } else {
        return(refdb)
    }
}

load_sqlite <- function(ah_id){
    path <- suppressMessages(ah[[ah_id]])
    conn <- dbConnect(SQLite(), path)
    return(conn)
}