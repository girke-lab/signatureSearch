#' @title gCMAP Search Method
#' @description 
#' It uses query signature to search against the reference database defined in 
#' the \code{\link{qSig}} object by gCMAP method, which is adapted from the 
#' gCMAP package (Sandmann et al., 2014)
#' @details 
#' The Bioconductor gCMAP (Sandmann et al. 2014) package provides access to a 
#' related but not identical implementation of the original CMAP algorithm 
#' proposed by Lamb et al. (2006). It uses as query a rank transformed GES and 
#' the reference database is composed of the labels of up and down regulated 
#' DEG sets. This is the opposite situation of the CMAP method, where the query 
#' is composed of the labels of up and down regulated DEGs and the database 
#' contains rank transformed GESs.
#' 
#' @section Column description:
#' Description of the score columns in the result table specific for gCMAP 
#' method:
#' \itemize{
#'     \item effect: Scaled bi-directional enrichment score corresponding to 
#'     the scaled_score under the CMAP result.
#'     \item nSet: Number of genes in the reference gene sets after applying
#'     the higher and lower cutoff.
#'     \item nFound: Number of genes in the reference gene sets that are 
#'     present in the query signature.
#'     \item signed: Whether the gene sets in the reference database have signs, 
#'     e.g. representing up and down regulated genes when computing scores.
#' }
#' Description of the other columns are available at the 'result' slot of the
#' \code{\link{gessResult}} object.
#' 
#' @param qSig \code{\link{qSig}} object defining the query signature, the GESS
#' method (should be 'gCMAP') and the path to the reference database.
#' @param higher The 'higher' threshold. If not 'NULL', genes with a score 
#' larger than 'higher' will be included in the gene set with sign +1. 
#' At least one of 'lower' and 'higher' must be specified.
#' @param lower The lower threshold. If not 'NULL', genes with a score smaller 
#' than 'lower' will be included in the gene set with sign -1. 
#' At least one of 'lower' and 'higher' must be specified.
#' @param chunk_size size of chunk per processing
#' @return \code{\link{gessResult}} object, the result table contains the 
#' search results for each perturbagen in the reference database ranked by 
#' their signature similarity to the query.
#' @seealso \code{\link{qSig}}, \code{\link{gessResult}}, \code{\link{gess}}
#' @references 
#' Sandmann, T., Kummerfeld, S. K., Gentleman, R., & Bourgon, R. 
#' (2014). gCMAP: user-friendly connectivity mapping with R. Bioinformatics , 
#' 30(1), 127–128. \url{https://doi.org/10.1093/bioinformatics/btt592}
#' @examples 
#' db_path <- system.file("extdata", "sample_db.h5", 
#'                        package = "signatureSearch")
#' library(signatureSearchData)
#' sample_db <- readHDF5chunk(db_path, colindex=1:100)
#' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
#' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
#' qsig_gcmap <- qSig(query=query_mat, gess_method="gCMAP", refdb=db_path)
#' gcmap <- gess_gcmap(qsig_gcmap, higher=1, lower=-1)
#' result(gcmap)
#' @export
#' 
gess_gcmap <- function(qSig, higher, lower, chunk_size=5000){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  #stopifnot(validObject(qSig))
  if(qSig@gess_method != "gCMAP"){
    stop(paste("The 'gess_method' slot of 'qSig' should be 'gCMAP'",
               "if using 'gess_gcmap' function"))
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
    cmap <- gCMAP::induceCMAPCollection(mat, higher=higher, lower=lower)
    c <- connectivity_score_raw(experiment=as.matrix(query), query=cmap)
    resultDF <- rbind(resultDF, data.frame(c))
  }
  ## Apply scaling of scores to full data set
  resultDF[,"effect"] <-.connnectivity_scale(resultDF$effect)
  resultDF <- resultDF[order(abs(resultDF$effect), decreasing=TRUE), ]
  row.names(resultDF) <- NULL
  resultDF <- sep_pcf(resultDF)
  # add target column
  target <- suppressMessages(get_targets(resultDF$pert))
  res <- left_join(resultDF, target, by=c("pert"="drug_name"))
  
  x <- gessResult(result = as_tibble(res),
                  query = qSig@query,
                  gess_method = qSig@gess_method,
                  refdb = qSig@refdb)
  return(x)
}

.connnectivity_scale <- function(scores) {
  p <- max(scores)
  q <- min(scores)
  ifelse(scores == 0, 0, ifelse( scores > 0, scores / p, -scores / q ))
}