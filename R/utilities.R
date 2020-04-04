sep_pcf <- function(res){
    new <- as.data.frame(t(vapply(seq_len(nrow(res)), function(i)
        unlist(strsplit(as.character(res$set[i]), "__")), 
        FUN.VALUE=character(3))), stringsAsFactors=FALSE)
    colnames(new) <- c("pert", "cell", "type")
    res <- cbind(new, res[,-1])
    return(res)
}

# readHDF5mat <- function(h5file, colindex=seq_len(10)) {
#     m <- h5read(h5file, "assay", index=list(NULL, colindex))
#     mycol <- h5read(h5file, "colnames", index=list(colindex, 1))
#     myrow <- h5read(h5file, "rownames")
#     h5closeAll()
#     rownames(m) <- as.character(myrow[,1])
#     colnames(m) <- as.character(mycol[,1])
#     return(m)
# }

# getH5dim <- function(h5file){
#     mat_dim <- h5ls(h5file)$dim[1]
#     mat_ncol <- as.numeric(gsub(".* x ","", mat_dim))
#     mat_nrow <- as.numeric(gsub(" x .*","", mat_dim))
#     return(c(mat_nrow, mat_ncol))
# }

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
        if(overwrite){
            create_empty_h5(h5file, delete_existing=TRUE)
        }
    } else {
        create_empty_h5(h5file, delete_existing=FALSE)
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
#' create_empty_h5(tmp_file, level=6)
#' @export
create_empty_h5 <- function(h5file, delete_existing=FALSE, level=6) {
    if(delete_existing) unlink(h5file)
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
#' create_empty_h5(tmp_file)
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
    if(printstatus) h5ls(h5file, all=TRUE)[c("dim", "maxdim")]
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
#' @param by_ncol number of columns to import in each iteration to limit 
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
#'         h5file=h5file, overwrite=TRUE)
#' @export
#' 
gctx2h5 <- function(gctx, cid, new_cid=cid, h5file, by_ncol=5000, 
                    overwrite=TRUE){
    cid_list <- suppressWarnings(
        split(cid, rep(seq_len(ceiling(length(cid)/by_ncol)), 
                       each=by_ncol)))
    new_cid_list <- suppressWarnings(
        split(new_cid, rep(seq_len(ceiling(length(new_cid)/by_ncol)), 
                           each=by_ncol)))
    if(file.exists(h5file)){
        if(overwrite){
            create_empty_h5(h5file, delete_existing=TRUE)
        }
    } else {
        create_empty_h5(h5file, delete_existing=FALSE)
    }
    lapply(seq_along(cid_list), function(i){
        mat <- parse_gctx(gctx, cid=cid_list[[i]], matrix_only=TRUE)
        mat <- mat@mat
        colnames(mat) <- new_cid_list[[i]]
        append2H5(x=mat, h5file, printstatus=FALSE)
    })    
    h5ls(h5file)
}

# readHDF5chunk <- function(h5file, colindex=seq_len(10), colnames=NULL) {
#     if(! is.null(colnames)){
#         all_trts <- h5read(h5file, "colnames", drop=TRUE)
#         colindex2 <- which(all_trts %in% colnames)
#         m <- h5read(h5file, "assay", index=list(NULL, colindex2))
#         colindex <- colindex2
#     } else {
#         m <- h5read(h5file, "assay", index=list(NULL, colindex))
#     }
#     mycol <- h5read(h5file, "colnames", index=list(colindex, 1))
#     myrow <- h5read(h5file, "rownames")
#     rownames(m) <- as.character(myrow[,1])
#     colnames(m) <- as.character(mycol[,1])
#     if(! is.null(colnames)){
#         m <- m[,colnames, drop=FALSE]
#     }
#     se <- SummarizedExperiment(assays=list(score=m))
#     h5closeAll()
#     return(se)
# }

#' @importFrom ExperimentHub ExperimentHub

determine_refdb <- function(refdb){
    eh <- suppressMessages(ExperimentHub())
    if(refdb=="cmap"){
        return(eh[["EH3223"]])
    }
    if(refdb=="cmap_expr"){
        return(eh[["EH3224"]])
    }
    if(refdb=="lincs"){
        return(eh[["EH3226"]])
    }
    if(refdb=="lincs_expr"){
        return(eh[["EH3227"]])
    } else {
        return(refdb)
    }
}

load_sqlite <- function(eh_id){
    eh <- suppressMessages(ExperimentHub())
    path <- suppressMessages(eh[[eh_id]])
    conn <- dbConnect(SQLite(), path)
    return(conn)
}

select_ont <- function(res, ont, GO_DATA){
    # Add and select ontology in res
    res <- add_GO_Ontology(res, GO_DATA)
    tmp_df <- result(res)
    colnames(tmp_df)[1] = "ont"
    tmp_df$ont = as.character(tmp_df$ont)
    rst(res) <- tmp_df
    if(ont != "ALL")
        rst(res) <- as_tibble(res[res$ont == ont, ])
    return(res)
}

#' Functionalities used to draw from reference database 
#' (e.g. lincs, lincs_expr) GESs of compound treatment(s) in cell types.
#' 
#' The GES could be genome-wide differential expression profiles (e.g. log2 
#' fold changes or z-scores) or normalized gene expression intensity values 
#' depending on the data type of \code{refdb} or n top up/down regulated DEGs
#' @title Drawn Query GES from Reference Database 
#' @rdname getSig
#' @param cmp character vector representing a list of compound name available 
#' in \code{refdb} for \code{getSig} function, or character(1) indicating a
#' compound name (e.g. vorinostat) for other functions
#' @param cell character(1) or character vector of the same length as cmp 
#' argument. It indicates cell type that the compound treated in
#' @param refdb character(1), one of "lincs", "lincs_expr", "cmap", "cmap_expr",
#' or path to the HDF5 file built from \code{\link{build_custom_db}} function
#' @return matrix representing genome-wide GES of the query compound(s) in cell
#' @examples 
#' refdb <- system.file("extdata", "sample_db.h5", package = "signatureSearch")
#' vor_sig <- getSig("vorinostat", "SKB", refdb=refdb)
#' @export
#' 
getSig <- function(cmp, cell, refdb){
    if(is.character(refdb)){
        refdb <- determine_refdb(refdb)
        refse <- SummarizedExperiment(HDF5Array(refdb, name="assay"))
        rownames(refse) <- HDF5Array(refdb, name="rownames")
        colnames(refse) <- HDF5Array(refdb, name="colnames")
        trt <- paste(cmp, cell, "trt_cp", sep="__")
        trt2 <- intersect(trt, colnames(refse))
        notin <- setdiff(trt, colnames(refse))
        if(length(notin) > 0){
            warning(length(notin), "/", length(trt), 
                    " teatments are not contained in refdb, ", 
                    "they are ignored!")
        }
        cmp_mat <- as.matrix(assay(refse[,trt2]))
        # sort decreasingly
        #cmp_mat2 <- apply(cmp_mat, 2, sort, decreasing=TRUE)
        return(cmp_mat)
    } else {
        message(paste("Please set refdb as one of",  
                "'lincs', 'lincs_expr', 'cmap' or 'cmap_expr', "),
                "or path to an HDF5 file representing reference database!")
    }
}

#' @rdname getSig
#' @param Nup integer(1). Number of most up-regulated genes to be subsetted
#' @param Ndown integer(1). Number of most down-regulated genes to be subsetted 
#' @return a list of up- and down-regulated gene label sets 
#' @examples 
#' vor_degsig <- getDEGSig("vorinostat", "SKB", refdb=refdb)
#' @export
getDEGSig <- function(cmp, cell, Nup=150, Ndown=150, refdb="lincs"){
    deprof <- suppressMessages(getSig(cmp, cell, refdb))
    deprof_sort <- apply(deprof, 2, sort, decreasing=TRUE)
    upset <- head(rownames(deprof_sort), Nup)
    downset <- tail(rownames(deprof_sort), Ndown)
    return(list(upset=upset, downset=downset))
}

#' @rdname getSig
#' @return a numeric matrix with one column representing gene expression values
#' drawn from \code{lincs_expr} db of the most up- and down-regulated genes. 
#' The genes were subsetted according to z-scores drawn from \code{lincs} db. 
#' @examples 
#' all_expr <- as.matrix(runif(1000, 0, 10), ncol=1)
#' rownames(all_expr) <- paste0('g', sprintf("%04d", 1:1000))
#' colnames(all_expr) <- "drug__cell__trt_cp"
#' de_prof <- as.matrix(rnorm(1000, 0, 3), ncol=1)
#' rownames(de_prof) <- paste0('g', sprintf("%04d", 1:1000))
#' colnames(de_prof) <- "drug__cell__trt_cp"
#' ## getSPsubSig internally uses deprof2subexpr function
#' ## sub_expr <- deprof2subexpr(all_expr, de_prof, Nup=150, Ndown=150)
#' @export
getSPsubSig <- function(cmp, cell, Nup=150, Ndown=150){
    query_expr <- suppressMessages(getSig(cmp, cell, refdb="lincs_expr"))
    query_prof <- suppressMessages(getSig(cmp, cell, refdb="lincs"))
    sub_expr <- deprof2subexpr(query_expr, query_prof, Nup=Nup, Ndown=Ndown)
    return(sub_expr)
}

deprof2subexpr <- function(all_expr, de_prof, Nup=150, Ndown=150){
    de_prof_sort <- apply(de_prof, 2, sort, decreasing=TRUE)
    upset <- head(rownames(de_prof_sort), Nup)
    downset <- tail(rownames(de_prof_sort), Ndown)
    sub_expr <- all_expr[c(upset, downset), , drop=FALSE]
    return(sub_expr)
}
    
trts_check <- function(ref_trts, full_trts){
    trts_valid <- intersect(ref_trts, full_trts)
    inval <- setdiff(ref_trts, full_trts)
    if(length(inval)>0){
        message(length(inval), 
        " treatments in ref_trts are not available in refdb, they are ignored!")
    }
    if(length(trts_valid)==0){
        stop("No ref_trts are available in refdb, ", 
             "the refdb is subsetted to empty!")
    }
    return(trts_valid)
}

#' Reduce number of characters for each element of a character vector by 
#' replacting the part that beyond Nchar (e.g. 50) character to '...'.  
#' @title Reduce Number of Character 
#' @param vec character vector to be reduced
#' @param Nchar integer, for each element in the vec, number of characters to 
#' remain
#' @return character vector after reducing
#' @examples 
#' vec <- c(strrep('a', 60), strrep('b', 30))
#' vec2 <- vec_char_redu(vec, Nchar=50)
#' @export
vec_char_redu <- function(vec, Nchar=50){
    vec <- as.character(vec)
    res <- vapply(vec, function(s){
        if(nchar(s) > Nchar){
            s2 <- substr(s, 1, 50)
            paste0(s2, "...")
        } else s
    }, FUN.VALUE=character(1))
    return(res)
}
