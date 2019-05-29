#' @title Correlation based method for GESS
#' @description 
#' It uses query signature to search against the reference database in the 
#' \code{\link{qSig}} object by measuring correlation coefficient.
#' @details 
#' The correlation coefficients can be used as a GESS method by searching 
#' with a query expression profile a database of expression profiles. 
#' Correlation-based queries were performed with genome-wide GEPs as well as 
#' with GEPs subset to the same query genes used for the set enrichment 
#' methods (\code{CMAP}, \code{LINCS} and \code{Fisher}). 
#' The latter situation makes the correlation-based results more comparable 
#' to the set enrichment methods by providing to each method a more equal 
#' amount of information than this is the case for the correlation method with 
#' genome-wide GEPs. 
#' 
#' Description of the score columns in the gess_cor tibble result:
#' \itemize{
#'     \item cor_score: correlation coefficiency.
#' }
#' @param qSig `qSig` object, The 'gess_method' slot should be 'Cor'. 
#' The reference database in the \code{qsig} could either store gene expression 
#' values or differential expression scores.
#' @param method One of 'spearman' (default), 'kendall', or 'pearson',
#' indicating which correlation coefficient to be used.
#' @param chunk_size size of chunk per processing
#' @return gessResult object, containing drugs in the reference database
#' ranked by their similarity to the query signature
#' @seealso \code{\link{qSig}}, \code{\link{gessResult}}, \code{\link{gess}}
#' @examples 
#' db_path <- system.file("extdata", "sample_db.h5", 
#'                        package = "signatureSearch")
#' library(signatureSearchData)
#' sample_db <- readHDF5chunk(db_path, colindex=1:100)
#' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
#' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
#' qsig_sp <- qSig(query = query_mat, gess_method = "Cor", refdb = db_path)
#' sp <- gess_cor(qSig=qsig_sp, method="spearman")
#' result(sp)
#' @export
gess_cor <- function(qSig, method, chunk_size=5000){
    if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
    #stopifnot(validObject(qSig))
    if(qSig@gess_method != "Cor"){
        stop("The 'gess_method' slot of 'qSig' should be 'Cor' 
             if using 'gess_cor' function")
  }
  query <- qSig@query
  db_path <- determine_refdb(qSig@refdb)
  mat_dim <- getH5dim(db_path)
  mat_ncol <- mat_dim[2]
  ceil <- ceiling(mat_ncol/chunk_size)
  resultDF <- data.frame()
  for(i in seq_len(ceil)){
    mat <- readHDF5mat(db_path,
                    colindex=(chunk_size*(i-1)+1):min(chunk_size*i, mat_ncol))
    cor_res <- cor_sig_search(query=query, refdb=mat, method=method)
    resultDF <- rbind(resultDF, data.frame(cor_res))
  }
  resultDF <- sep_pcf(resultDF)
  resultDF <- resultDF[order(abs(resultDF$cor_score), decreasing = TRUE), ]
  row.names(resultDF) <- NULL
  # add target column
  target <- suppressMessages(get_targets(resultDF$pert))
  res <- left_join(resultDF, target, by=c("pert"="drug_name"))
  
  x <- gessResult(result = as_tibble(res),
                  query = qSig@query,
                  gess_method = qSig@gess_method,
                  refdb = qSig@refdb)
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


