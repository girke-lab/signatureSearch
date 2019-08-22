#' @title Correlation-based Search Method
#' @description 
#' Correlation-based similarity metrics, such as Spearman or Pearson coefficients,
#' can be used as Gene Expression Signature Search (GESS) methods. As non-set-based 
#' methods, they require quantitative gene expression values for both the query and
#' the database entries, such as normalized intensities or read counts from
#' microarrays or RNA-Seq experiments, respectively.
#' @details 
#' For correlation searches to work, it is important that both the query and
#' reference database contain the same type of gene identifiers. The expected data
#' structure of the query is a matrix with a single numeric column and the gene
#' labels (e.g. Entrez Gene IDs) in the row name slot. For convenience, the
#' correlation-based searches can either be performed with the full set of genes
#' represented in the database or a subset of them. The latter can be useful to
#' focus the computation for the correlation values on certain genes of interest
#' such as a DEG set or the genes in a pathway of interest. For comparing the
#' performance of different GESS methods, it can also be advantageous to subset
#' the genes used for a correlation-based search to same set used in a set-based
#' search, such as the up/down DEGs used in a LINCS GESS. This way the search
#' results of correlation- and set-based methods can be more comparable because
#' both are provided with equivalent information content.
#' 
#' @section Column description:
#' Descriptions of the columns specific to the corrleation-based GESS method are
#' given below. Note, the additional columns, those that are common among the 
#' GESS methods, are described in the help file of the \code{gessResult} object.
#' \itemize{
#'     \item cor_score: Correlation coefficient based on the method defined in 
#'     the \code{gess_cor} function.
#' }
#' 
#' @param qSig \code{\link{qSig}} object defining the query signature including
#' the GESS method (should be 'Cor') and the path to the reference database. For
#' details see help of \code{qSig} and \code{qSig-class}.
#' @param method One of 'spearman' (default), 'kendall', or 'pearson',
#' indicating which correlation coefficient to use.
#' @param chunk_size number of database entries to process per iteration to limit
#' memory usage of search.
#' @return \code{\link{gessResult}} object, the result table contains the 
#' search results for each perturbagen in the reference database ranked by 
#' their signature similarity to the query.
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


