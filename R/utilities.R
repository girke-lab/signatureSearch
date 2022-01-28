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
#' @param name The name of the dataset in the HDF5 file. The default is write
#' the score matrix (e.g. z-score, logFC) to the 'assay' dataset, users could 
#' also write the adjusted p-value or FDR matrix to the 'padj' dataset by 
#' setting the \code{name} as 'padj'.
#' @param overwrite TRUE or FALSE, whether to overwrite or append
#' matrix to an existing 'h5file'
#' @return HDF5 file containing exported matrix
#' @examples 
#' mat <- matrix(rnorm(12), nrow=3, dimnames=list(
#'               paste0("r",1:3), paste0("c",1:4)))
#' h5file <- tempfile(fileext=".h5")
#' matrix2h5(matrix=mat, h5file=h5file, overwrite=TRUE)
#' @export
matrix2h5 <- function(matrix, h5file, name="assay", overwrite=TRUE){
    create_empty_h5(h5file, delete_existing=overwrite)
    append2H5(matrix, h5file, name=name)
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
    if(! file.exists(h5file)){
        h5createFile(file=h5file)
        h5createDataset(h5file, "assay", c(0,0), c(H5Sunlimited(), H5Sunlimited()), 
                        chunk=c(12328,1), level=level)
        h5createDataset(h5file, "padj", c(0,0), c(H5Sunlimited(), H5Sunlimited()), 
                        chunk=c(12328,1), level=level)
        h5createDataset(h5file, "colnames", c(0,1), c(H5Sunlimited(), 1), 
                        storage.mode='character', size=200, level=level)
        h5createDataset(h5file, "rownames", c(0,1), c(H5Sunlimited(), 1), 
                        storage.mode='character', size=40, level=level)
    }
    h5closeAll()
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
#' @param name The name of the dataset in the HDF5 file.
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
append2H5 <- function(x, h5file, name="assay", printstatus=TRUE) {
    status <- h5ls(h5file)[c("name", "dim")]
    rowstatus <- as.numeric(gsub(" x \\d{1,}$", "", 
                                 status[status$name==name, "dim"]))
    colstatus <- as.numeric(gsub("^\\d{1,} x ", "", 
                                 status[status$name==name, "dim"]))
    nrows <- nrow(x) 
    ncols <- colstatus + ncol(x)
    h5set_extent(h5file, name, c(nrows, ncols))
    h5write(x, h5file, name, index=list(seq_len(nrows), (colstatus+1):ncols))
    # only stores rownames and colnames of 'assay' dataset
    if(name=="assay"){
        h5set_extent(h5file, "colnames", c(ncols,1))
        h5write(colnames(x), h5file, "colnames", index=list((colstatus+1):ncols, 1))
        if(any(duplicated(h5read(h5file, "colnames")[,1]))) 
            warning("Column names contain duplicates!")
        
        h5set_extent(h5file, "rownames", c(nrows,1))
        h5write(rownames(x), h5file, "rownames", index=list(seq_len(nrows), 1))
        if(any(duplicated(h5read(h5file, "rownames")[,1]))) 
            warning("Row names contain duplicates!")
    }
    
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
    create_empty_h5(h5file, delete_existing=overwrite)
    lapply(seq_along(cid_list), function(i){
        mat <- parse_gctx(gctx, cid=cid_list[[i]], matrix_only=TRUE)
        mat <- mat@mat
        colnames(mat) <- new_cid_list[[i]]
        append2H5(x=mat, h5file, printstatus=FALSE)
    })    
    h5ls(h5file)
}

#' Calculate Mean Expression Values of LINCS Level 3 Data
#' 
#' Function calculates mean expression values for replicated samples of LINCS
#' Level 3 data that have been treated by the same compound in the same cell type
#' at a chosen concentration and treatment time. Usually, the function is used
#' after filtering the Level 3 data with \code{inst_filter}. The results (here
#' matrix with mean expression values) are saved to an HDF5 file. The latter is
#' referred to as the `lincs_expr` database.
#' @param gctx character(1), path to the LINCS Level 3 gctx file
#' @param inst tibble, LINCS Level 3 instances after filtering for specific concentrations and times
#' @param h5file character(1), path to the destination HDF5 file
#' @param chunksize number of columns of the matrix to be processed at a time to limit memory usage
#' @param overwrite TRUE or FALSE, whether to overwrite or append data to an
#' existing 'h5file'
#' @return HDF5 file, representing the \code{lincs_expr} database
#' @examples
#' gctx <- system.file("extdata", "test_sample_n2x12328.gctx", package="signatureSearch")
#' h5file <- tempfile(fileext=".h5")
#' inst <- data.frame(inst_id=c("ASG001_MCF7_24H:BRD-A79768653-001-01-3:10",
#'                              "CPC012_SKB_24H:BRD-K81418486:10"), 
#'     pert_cell_factor=c('sirolimus__MCF7__trt_cp', 'vorinostat__SKB__trt_cp'))
#' meanExpr2h5(gctx, inst, h5file, overwrite=TRUE)
#' @export
#' 
meanExpr2h5 <- function(gctx, inst, h5file, chunksize=2000, overwrite=TRUE){
    ## Get mean expression values of replicate samples from the same drug treatment in cells
    pert_cell_list <- split(inst$inst_id, inst$pert_cell_factor)
    chunk_list <- suppressWarnings(
        split(pert_cell_list, rep(seq_len(ceiling(length(pert_cell_list)/chunksize)), each=chunksize)))
    ## Creat an empty h5file
    if(file.exists(h5file)){
        if(overwrite){
            create_empty_h5(h5file, delete_existing=TRUE)
        }
    } else {
        create_empty_h5(h5file, delete_existing=FALSE)
    }
    ## append mean expr mat in each chunk to h5file
    lapply(chunk_list, function(pcl){
        cid_all <- unlist(pcl)
        mat_all <- parse_gctx(gctx, cid=cid_all, matrix_only=TRUE)
        mat_all <- mat_all@mat
        expr <- vapply(pcl, function(cid, mat_all){
            rowMeans(as.matrix(mat_all[,cid]))
        }, mat_all=mat_all, FUN.VALUE=numeric(12328))
        append2H5(x=expr, h5file, printstatus=FALSE)
    })
    h5ls(h5file)
}

#' Read gene sets from large gmt file in batches, convert the gene sets to
#' 01 matrix and write the result to an HDF5 file.
#' @title Convert GMT to HDF5 File
#' @param gmtfile character(1), path to gmt file containing gene sets
#' @param dest_h5 character(1), path of the hdf5 destination file
#' @param by_nset number of gene sets to import in each iteration to limit 
#' memory usage
#' @param overwrite TRUE or FALSE, whether to overwrite or to append to 
#' existing 'h5file'
#' @return HDF5 file
#' @examples 
#' gmt <- system.file("extdata", "test_gene_sets_n4.gmt", 
#'         package="signatureSearch")
#' h5file <- tempfile(fileext=".h5")
#' gmt2h5(gmtfile=gmt, dest_h5=h5file, overwrite=TRUE)
#' @export
gmt2h5 <- function(gmtfile, dest_h5, by_nset=5000, overwrite=FALSE){
    # get number of lines of gmtfile
    wc <- system(paste("wc -l", gmtfile), intern = TRUE)
    nline <- as.numeric(gsub(pattern=gmtfile, "", wc,fixed = TRUE))
    ceil <- ceiling(nline/by_nset)
    
    if(file.exists(dest_h5) & !overwrite){
        message(paste(dest_h5, "already exists!"))
    } else {
        create_empty_h5(dest_h5, delete_existing=TRUE)
    }
    
    # get all gene identifiers in gmtfile
    all_genes <- NULL
    for(i in seq_len(ceil)){
        gene_sets <- suppressWarnings(read_gmt(gmtfile, 
                            start=by_nset*(i-1)+1, end=by_nset*i))
        tmp <- unique(unlist(gene_sets))
        all_genes <- unique(c(all_genes, tmp))
    }
            
    # read in gene sets in batches
    for(i in seq_len(ceil)){
        gene_sets <- suppressWarnings(read_gmt(gmtfile, 
                        start=by_nset*(i-1)+1, end=by_nset*i))
        # transform gene sets to sparseMatrix
        mat <- gs2mat(gene_sets)
        miss_genes <- setdiff(all_genes, rownames(mat))
        patch <- matrix(0, nrow=length(miss_genes), ncol=ncol(mat))
        rownames(patch) <- miss_genes
        colnames(patch) <- colnames(mat)
        mat <- rbind(mat, patch)
        append2H5(x=mat[all_genes,], dest_h5, printstatus=FALSE)
    }
    h5ls(dest_h5)
}

gs2mat <- function(gene_sets){
    gsc <- GeneSetCollection(mapply(function(geneIds, setId) {
        GeneSet(geneIds, geneIdType=EntrezIdentifier(),
                setName=setId)
    }, gene_sets, names(gene_sets)))
    mat <- incidence(gsc)
    mat <- Matrix::t(mat)
    return(mat)
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

select_ont <- function(res, ont, GO_DATA){
    # Add and select ontology in res
    res <- add_GO_Ontology(res, GO_DATA)
    tmp_df <- result(res)
    colnames(tmp_df)[1] = "ont"
    tmp_df$ont = as.character(tmp_df$ont)
    result(res) <- tmp_df
    if(ont != "ALL")
        result(res) <- as_tibble(res[res$ont == ont, ])
    return(res)
}

#' Functionalities used to draw from reference database 
#' (e.g. lincs, lincs_expr) GESs of compound treatment(s) in cell types.
#' 
#' The GES could be genome-wide differential expression profiles (e.g. log2 
#' fold changes or z-scores) or normalized gene expression intensity values 
#' depending on the data type of \code{refdb} or n top up/down regulated DEGs
#' @title Draw GESs from Reference Database 
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
#' @param higher numeric(1), the upper threshold of defining DEGs. 
#' At least one of 'lower' and 'higher' must be specified.
#' If \code{Nup} or \code{Ndown} arguments are defined, it will be ignored. 
#' @param lower numeric(1), the lower threshold of defining DEGs. 
#' At least one of 'lower' and 'higher' must be specified.
#' If \code{Nup} or \code{Ndown} arguments are defined, it will be ignored. 
#' @param padj numeric(1), cutoff of adjusted p-value or false discovery rate (FDR)
#' of defining DEGs if the reference HDF5 database contains the p-value matrix 
#' stored in the dataset named as 'padj'.
#' If \code{Nup} or \code{Ndown} arguments are defined, it will be ignored. 
#' @return a list of up- and down-regulated gene label sets 
#' @examples 
#' vor_degsig <- getDEGSig(cmp="vorinostat", cell="SKB", Nup=150, Ndown=150,
#'                         refdb=refdb)
#' @export
getDEGSig <- function(cmp, cell, Nup=NULL, Ndown=NULL, 
                      higher=NULL, lower=NULL, padj=NULL, refdb="lincs"){
    deprof <- suppressMessages(getSig(cmp, cell, refdb))
    deprof_sort <- apply(deprof, 2, sort, decreasing=TRUE)
    if(!is.null(Nup) & is.null(Ndown)){
        Ndown = 0
    }
    if(is.null(Nup) & !is.null(Ndown)){
        Nup = 0
    }
    if(! is.null(Nup) & ! is.null(Ndown)){
        upset <- head(rownames(deprof_sort), Nup)
        downset <- tail(rownames(deprof_sort), Ndown)
        return(list(upset=upset, downset=downset))
    } 
    if(!is.null(higher) & is.null(lower)){
        up <- rownames(deprof_sort[deprof_sort >= higher,,drop=FALSE])
        if(!is.null(padj)){
            pdegs <- getPCut(cmp, cell, refdb, padj)
            up <- up[up %in% pdegs]
        }
        return(list(upset=up, downset=character()))
    }
    if(is.null(higher) & !is.null(lower)){
        down <- rownames(deprof_sort[deprof_sort <= lower,,drop=FALSE])
        if(!is.null(padj)){
            pdegs <- getPCut(cmp, cell, refdb, padj)
            down <- down[down %in% pdegs]
        }
        return(list(upset=character(), downset=down))
    }
    if(!is.null(higher) & !is.null(lower)){
        up <- rownames(deprof_sort[deprof_sort >= higher,,drop=FALSE])
        down <- rownames(deprof_sort[deprof_sort <= lower,,drop=FALSE])
        if(!is.null(padj)){
            pdegs <- getPCut(cmp, cell, refdb, padj)
            up <- up[up %in% pdegs]
            down <- down[down %in% pdegs]
        }
        return(list(upset=up, downset=down))
    }
    stop("You need to either set 'Nup', 'Ndown' as integers or 
    set 'higher', 'lower', 'padj' as numeric values")
}

getPCut <- function(cmp, cell, refdb, padj){
    db_path <- determine_refdb(refdb)
    pmat <- HDF5Array(db_path, name="padj")
    rownames(pmat) <- as.character(HDF5Array(db_path, name="rownames"))
    colnames(pmat) <- as.character(HDF5Array(db_path, name="colnames"))
    pval <- as.matrix(pmat[, paste(cmp, cell, "trt_cp", sep="__"), drop=FALSE])
    degs <- rownames(pval[pval <= padj,,drop=FALSE])
    return(degs)
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

#' Show Reduced Targets
#' 
#' Reduce number of targets for each element of a character vector by 
#' replacting the targets that beyond Ntar to '...'.  
#' @param vec character vector, each element composed by a list of targets 
#' symbols separated by '; '
#' @param Ntar integer, for each element in the vec, number of targets to show
#' @return character vector after reducing
#' @examples 
#' vec <- c("t1; t2; t3; t4; t5; t6", "t7; t8")
#' vec2 <- tarReduce(vec, Ntar=5)
#' @export
#' 
tarReduce <- function(vec, Ntar=5){
    singleTarShot <- function(s, Ntar){
        temp <- strsplit(s, ';\\s?')[[1]]
        if(length(temp)>Ntar){
            return(paste0(paste(temp[seq_len(Ntar)], collapse = "; "), "; ..."))
        } else {
            return(s)
        }
    }
    return(sapply(vec, singleTarShot, Ntar, USE.NAMES = FALSE))
}

load_OrgDb <- function(OrgDb){
    if(is(OrgDb, "character")){
        if(! require(OrgDb, character.only = TRUE)){
            stop(paste(OrgDb, "package need to be installed to use this function"))
        }
        OrgDb <- eval(parse(text = OrgDb))
    }
    return(OrgDb)
}

#' Read in gene set information from .gmt files
#' 
#' This function reads in and parses information from the MSigDB's .gmt files. 
#' Pathway information will be returned as a list of gene sets.
#' 
#' The .gmt format is a tab-delimited list of gene sets, where each line is a 
#' separate gene set. The first column must specify the name of the gene set, 
#' and the second column is used for a short description (which this function 
#' discards). For complete details on the .gmt format, refer to the Broad 
#' Institute's Data Format's page 
#' \url{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}.
#' 
#' @param file The .gmt file to be read
#' @param start integer(1), read the gmt file from start line
#' @param end integer(1), read the gmt file to the end line, the default -1 
#' means read to the end 
#' @return A list, where each index represents a separate gene set.
#' @section Warning:
#' The function does not check that the file is correctly formatted, and may 
#' return incorrect or partial gene sets, e.g. if the first two columns are 
#' omitted. Please make sure that files are correctly formatted before reading 
#' them in using this function.
#' @examples
#' gmt_path <- system.file("extdata/test_gene_sets_n4.gmt", package="signatureSearch")
#' geneSets <- read_gmt(gmt_path)
#' @export
#' 
read_gmt <- function(file, start=1, end=-1){
    if (!grepl("\\.gmt$", file)[1]) {
        stop("Pathway information must be a .gmt file")
    }
    geneSetDB = scan(file, what="character", n=end-start+1, skip=start-1, 
                     sep="\n", quiet=TRUE)
    geneSetDB = suppressWarnings(strsplit(geneSetDB, "\t"))
    names(geneSetDB) = sapply(geneSetDB, "[", 1)
    geneSetDB = lapply(geneSetDB, "[", -1:-2)
    geneSetDB = lapply(geneSetDB, function(x) {
        x[which(x != "")]
    })
    geneSetDB <- geneSetDB[sapply(geneSetDB, length)>0 & !is.na(names(geneSetDB))]
    return(geneSetDB)
}

##' Mapping 'itemID' column in the FEA enrichment result table from Entrez ID
##' to gene Symbol
##'
##' @title Set Readable
##' @param tb tibble object, enrichment result table
##' @param OrgDb character(1), 'org.Hs.eg.db' for human
##' @param keyType character(1), keyType of gene
##' @param geneCol character(1), name of the column in 'tb' containing gene
##' Entrez ids separated by '/' to be converted to gene Symbol
##' @return tibble Object
##' @importFrom stats setNames
##' @importFrom AnnotationDbi columns
##' @examples 
##' data(drugs10)
##' res <- tsea_dup_hyperG(drugs=drugs10, type="Reactome", pvalueCutoff=1, 
##'                        qvalueCutoff=1)
##' res_tb <- set_readable(result(res))
##' @export
set_readable <- function(tb, OrgDb="org.Hs.eg.db", keyType="ENTREZID", geneCol="itemID"){
    OrgDb <- load_OrgDb(OrgDb)
    if (!'SYMBOL' %in% columns(OrgDb)) {
        warning("Fail to convert input geneID to SYMBOL since no SYMBOL information available in the provided OrgDb...")
    }
    
    if(! geneCol %in% colnames(tb)){
        stop("The column name defined under 'geneCol' must be contained in tb!")
    }
    gc <- setNames(strsplit(tb[[geneCol]], "/", fixed=TRUE), tb$ID)
    genes <- unique(unlist(gc))
    gn <- EXTID2NAME(OrgDb, genes, keyType)
    gc <- lapply(gc, function(i) gn[i])
    geneID <- sapply(gc, paste0, collapse="/")
    
    # gc_sym <- sapply(gc, function(genes){
    #     gn_df <- suppressMessages(select(OrgDb, keys=genes, keytype=keyType, 
    #                                      columns="SYMBOL"))
    #     gn_sym <- paste0(na.omit(gn_df$SYMBOL), collapse="/")
    #     return(gn_sym)})
    
    tb[[geneCol]] <- geneID
    return(tb)
}

#' Add PCID to drug data frame
#' 
#' This function can be used to add the \code{PCIDss} (PubChem CID column added
#' from signatureSearch package) column to a data frame that have a column 
#' store compound names. The compound name to PubChem CID annotation is obtained
#' from \code{lincs_pert_info} in 2017.
#' 
#' @param df data frame or tibble object
#' @param drug_col name of the column that store compound names in df
#' @return tibble object with an added PCIDss column
#' @importFrom dplyr relocate
#' @examples 
#' data("lincs_pert_info")
#' # gess_tb2 <- add_pcid(gess_tb)
#' @export
add_pcid <- function(df, drug_col="pert"){
    data("lincs_pert_info", envir=environment())
    pert2 <- lincs_pert_info[, c("pert_iname", "pubchem_cid")]
    colnames(pert2) <- c("pert_iname", "PCIDss")
    join_cols = "pert_iname"
    names(join_cols) <- drug_col
    df %<>% left_join(pert2, by=join_cols)
    return(df)
}

#' Add MOA annotation to drug data frame
#' 
#' The MOA annotation is a list of MOA name to drug name mappings. 
#' This functions add the MOA column to data frame when data frame have a 
#' column with compound names
#' 
#' @param df data frame that must contains a column with compound names
#' @param drug_col character (1), name of the column that stores compound names
#' @param moa_list a list object of MOA name (e.g. HDAC inhibitor) to compound
#' name mappings
#' @return data frame with an added MOAss column  
#' @examples 
#' data("clue_moa_list")
#' df <- data.frame(pert=c("vorinostat", "sirolimus"), annot1=c("a", "b"), 
#'                  annot2=1:2)
#' addMOA(df, "pert", clue_moa_list) 
#' @export
addMOA <- function(df, drug_col, moa_list){
    drug2moa <- list_rev(moa_list)
    d2m_vec <- sapply(drug2moa, paste, collapse="; ")
    df$MOAss <- d2m_vec[df[[drug_col]]]
    return(df)
}

#' Add Compound Annotation Info to GESS Result Table
#' 
#' This function supports adding customized compound annotation table to the 
#' GESS result table if provided. It then automatically adds the target gene
#' symbols, MOAs and PubChem CID columns (\code{t_gn_sym}, \code{MOAss}, 
#' \code{PCIDss}) if the table contains a column that stores compound names.
#'  
#' @param gess_tb tibble or data.frame object of GESS result, 
#' can be accessed by the \code{result} method on the \code{gessResult} object 
#' from \code{gess_*} functions. Or a customized data frame that contains a 
#' \code{pert} column that stores compound id or name.
#' @param refdb character(1), reference database that can be accessed by 
#' the \code{refdb} method on the \code{gessResult} object. If \code{gess_tb}
#' is a customized table, \code{refdb} can be just set as 'custom'.
#' @param cmp_annot_tb data.frame or tibble of compound annotation table. 
#' This table contains annotation information for compounds stored under
#' \code{pert} column of \code{gess_tb}. Set to \code{NULL} if not available.
#' @param by character(1), column name in \code{cmp_annot_tb} that can be merged
#' with \code{pert} column in \code{gess_tb} if \code{refdb} is not set 
#' as 'lincs2', otherwise, it is the column name in \code{cmp_annot_tb} that can 
#' be merged with \code{pert_id} column in the GESS result table. 
#' If \code{cmp_annot_tb} is NULL, \code{by} is ignored.
#' @param cmp_name_col character(1), column name in \code{gess_tb} or 
#' \code{cmp_annot_tb} that store compound names. If there is no compound name 
#' column, set to \code{NULL}. If \code{cmp_name_col} is available, 
#' three additional columns (t_gn_sym, MOAss, PCIDss) are automatically added 
#' by using \code{\link{get_targets}}, CLUE touchstone compound MOA annotation, 
#' and 2017 \code{lincs_pert_info} annotation table, respectively as 
#' annotation sources. t_gn_sym: target gene symbol, MOAss: MOA annotated from 
#' signatureSearch, PCIDss: PubChem CID annotated from signatureSearch.
#' @return tibble of \code{gess_tb} with target, MOA, PubChem CID annotations
#' and also merged with user provided compound annotation table. 
#' @importFrom dplyr tibble
#' @examples 
#' gess_tb <- data.frame(pert=c("vorinostat", "sirolimus", "estradiol"),
#'                   cell=c("SKB", "SKB", "MCF7"),
#'                   NCS=runif(3))
#' cmp_annot_tb <- data.frame(pert_name=c("vorinostat", "sirolimus", "estradiol"),
#'                            prop1=c("a", "b", "c"),
#'                            prop2=1:3)
#' addGESSannot(gess_tb, "custom", cmp_annot_tb, by="pert_name", 
#'              cmp_name_col="pert")
#' @export
addGESSannot <- function(gess_tb, refdb, cmp_annot_tb=NULL, by="pert", 
                         cmp_name_col="pert"){
    names(by) <- "pert"
    lastcol <- colnames(gess_tb)[ncol(gess_tb)]
    if(refdb == "lincs2"){
        # rename 'pert' as 'pert_id'; add 'pert' column
        data("lincs_pert_info2", envir=environment())
        gess_tb %<>% left_join(lincs_pert_info2[,c("pert_id", "pert_iname")],
                               by=c("pert"="pert_id")) %>% 
            relocate("pert_iname", .after="pert") %>% 
            dplyr::rename(c("pert_id"="pert", "pert"="pert_iname"))
        names(by) <- "pert_id"
    }
    if(!is.null(cmp_annot_tb)){
        if(! by %in% colnames(cmp_annot_tb)){
            warning("'by' argument needs to be a column name that is in 'cmp_annot_tb' so that 
  this column can be merged with 'pert' column in GESS results. 
  Now, the compound annotation table will not be added to GESS result table.")
        } else {
            res <- left_join(tibble(gess_tb), cmp_annot_tb, by=by)
        }
    } else {
        res <- gess_tb
    }
    
    if(is.null(cmp_name_col)){
        return(tibble(res)) 
    } else if(! cmp_name_col %in% colnames(res)){
                cmp_name_col <- "pert"
    }
    
    # Add t_gn_sym
    target <- suppressMessages(get_targets(res[[cmp_name_col]]))
    join_cols <- "drug_name"; names(join_cols) <- cmp_name_col
    res %<>% left_join(target, by=join_cols)
    # Add MOAss
    data("clue_moa_list", envir=environment())
    res <- addMOA(res, cmp_name_col, clue_moa_list)
    # Add PCIDss
    res <- add_pcid(res, cmp_name_col)
    res %<>% relocate(t_gn_sym:PCIDss, .after=all_of(lastcol))
    return(tibble(res))
}

#' Named list to data frame 
#' 
#' Convert a list with names that have one to many mapping relationships to
#' a data.frame of two columns, one column is names, the other column is the
#' unlist elements
#' 
#' @param list input list with names slot
#' @param colnames character vector of length 2, indicating the column names
#' of the returned data.frame
#' @return data.frame
#' @examples 
#' list <- list("n1"=c("e1", "e2", "e4"), "n2"=c("e3", "e5"))
#' list2df(list, colnames=c("name", "element"))
#' @export
list2df <- function(list, colnames){
    df_list <- lapply(names(list), function(x){
        return(data.frame(x=x, y=list[[x]]))
    })
    df <- do.call(rbind, df_list)
    colnames(df) <- colnames
    return(unique(df))
}

#' Reverse list 
#' 
#' Reverse list from list names to elements mapping to elements to names mapping.
#' 
#' @param list input list with names slot
#' @return list
#' @examples 
#' list <- list("n1"=c("e1", "e2", "e4"), "n2"=c("e1", "e5"))
#' list_rev(list)
#' @export
list_rev <- function(list){
    df <- list2df(list, colnames=c("name", "ele"))
    listr <- split(df$name, df$ele)
    listr2 <- sapply(listr, unique, simplify=FALSE)
    return(listr2)
}

## get rid of "Undefined global functions or variables" note
globalVariables(c("PCID", "pubchem_cid", "lincs_pert_info",
                  "PCIDss", "all_of", "lincs_pert_info2", "t_gn_sym")) 
