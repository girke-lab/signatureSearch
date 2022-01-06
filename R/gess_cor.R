#' @rdname gess
#' @description 
#' Correlation-based similarity metrics, such as Spearman or Pearson 
#' coefficients, can be used as Gene Expression Signature Search (GESS) methods.
#' As non-set-based methods, they require quantitative gene expression values 
#' for both the query and the database entries, such as normalized intensities 
#' or read counts from microarrays or RNA-Seq experiments, respectively.
#' @details 
#' For correlation searches to work, it is important that both the query and
#' reference database contain the same type of gene identifiers. The expected 
#' data structure of the query is a matrix with a single numeric column and the 
#' gene labels (e.g. Entrez Gene IDs) in the row name slot. For convenience, the
#' correlation-based searches can either be performed with the full set of genes
#' represented in the database or a subset of them. The latter can be useful to
#' focus the computation for the correlation values on certain genes of interest
#' such as a DEG set or the genes in a pathway of interest. For comparing the
#' performance of different GESS methods, it can also be advantageous to subset
#' the genes used for a correlation-based search to same set used in a set-based
#' search, such as the up/down DEGs used in a LINCS GESS. This way the search
#' results of correlation- and set-based methods can be more comparable because
#' both are provided with equivalent information content.
#' @examples 
#' 
#' ######## Correlation-based GESS method #########
#' # qsig_sp <- qSig(query=query_mat, gess_method="Cor", refdb=db_path)
#' # sp <- gess_cor(qSig=qsig_sp, method="spearman")
#' # result(sp)
#' @export
gess_cor <- function(qSig, method="spearman", 
                     chunk_size=5000, ref_trts=NULL, workers=1, 
                     cmp_annot_tb=NULL, by="pert", cmp_name_col="pert"){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  #stopifnot(validObject(qSig))
  if(gm(qSig) != "Cor"){
      stop("The 'gess_method' slot of 'qSig' should be 'Cor' 
           if using 'gess_cor' function")
  }
  query <- qr(qSig)
  db_path <- determine_refdb(refdb(qSig))
  
  ## calculate cs_raw of query to blocks (e.g., 5000 columns) of full refdb
  full_mat <- HDF5Array(db_path, "assay")
  rownames(full_mat) <- as.character(HDF5Array(db_path, "rownames"))
  colnames(full_mat) <- as.character(HDF5Array(db_path, "colnames"))
  
  if(! is.null(ref_trts)){
      trts_valid <- trts_check(ref_trts, colnames(full_mat))
      full_mat <- full_mat[, trts_valid]
  }
  
  full_dim <- dim(full_mat)
  full_grid <- colAutoGrid(full_mat, ncol=min(chunk_size, ncol(full_mat)))
  ### The blocks in 'full_grid' are made of full columns 
  nblock <- length(full_grid) 
  resultDF <- bplapply(seq_len(nblock), function(b){
    ref_block <- read_block(full_mat, full_grid[[b]])
    cor_res <- cor_sig_search(query=query, refdb=ref_block, method=method)
    return(data.frame(cor_res))}, BPPARAM = MulticoreParam(workers = workers))
  resultDF <- do.call(rbind, resultDF)
  
  # mat_dim <- getH5dim(db_path)
  # mat_ncol <- mat_dim[2]
  # ceil <- ceiling(mat_ncol/chunk_size)
  # resultDF <- data.frame()
  # for(i in seq_len(ceil)){
  #   mat <- readHDF5mat(db_path,
  #                   colindex=(chunk_size*(i-1)+1):min(chunk_size*i, mat_ncol))
  #   cor_res <- cor_sig_search(query=query, refdb=mat, method=method)
  #   resultDF <- rbind(resultDF, data.frame(cor_res))
  # }

  resultDF <- resultDF[order(abs(resultDF$cor_score), decreasing = TRUE), ]
  row.names(resultDF) <- NULL
  res <- sep_pcf(resultDF)
  # add compound annotations
  res <- addGESSannot(res, refdb(qSig), cmp_annot_tb, by, cmp_name_col)
  x <- gessResult(result = res,
                  query = qr(qSig),
                  gess_method = gm(qSig),
                  refdb = refdb(qSig))
  return(x)
}

cor_sig_search <- function(query, refdb, method){
  res_list <- NULL
  # make sure rownames of query and refdb are the same
  common_gene <- intersect(rownames(query), rownames(refdb))
  query2 <- as.matrix(query[common_gene,])
  colnames(query2) <- colnames(query)
  refdb <- refdb[common_gene,]
  for(i in seq_len(ncol(query2))){
    cor <- as.numeric(cor(query2[,i], refdb, method = method))
    names(cor) <- colnames(refdb)
    trend=cor
    trend[cor>=0]="up"
    trend[cor<0]="down"
    res <- data.frame(set=names(cor), trend=trend, cor_score = cor, 
                      stringsAsFactors = FALSE)
    res_list <- c(res_list, list(res))
  }
  names(res_list) <- colnames(query)
  if(length(res_list)==1) return(res_list[[1]])
  return(res_list)
}


