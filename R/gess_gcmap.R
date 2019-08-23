#' @title gCMAP Search Method
#' @description 
#' Adapts the Gene Expression Signature Search (GESS) method from the gCMAP
#' package (Sandmann et al. 2014) to make it compatible with the database
#' containers and methods defined by \code{signatureSearch}. The specific GESS
#' method, called gCMAP, uses as query a rank transformed GES and the reference
#' database is composed of the labels of up and down regulated DEG sets.
#' @details 
#' The Bioconductor gCMAP (Sandmann et al. 2014) package provides access to a 
#' related but not identical implementation of the original CMAP algorithm 
#' proposed by Lamb et al. (2006). It uses as query a rank transformed GES and 
#' the reference database is composed of the labels of up and down regulated 
#' DEG sets. This is the opposite situation of the orignal CMAP method from 
#' Lamb et al (2006), where the query is composed of the labels of up and down 
#' regulated DEGs and the database contains rank transformed GESs.
#' 
#' @section Column description:
#' Descriptions of the columns specific to the gCMAP method are given below. 
#' Note, the additional columns, those that are common among the GESS methods, 
#' are described in the help file of the \code{gessResult} object.
#' 
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
#' 
#' @param qSig \code{\link{qSig}} object defining the query signature including
#' the GESS method (should be 'gCMAP') and the path to the reference database.
#' For details see help of \code{qSig} and \code{qSig-class}.
#' @param higher The 'upper' threshold. If not 'NULL', genes with a score 
#' larger than 'higher' will be included in the gene set with sign +1. 
#' At least one of 'lower' and 'higher' must be specified.
#' @param lower The lower threshold. If not 'NULL', genes with a score smaller 
#' than 'lower' will be included in the gene set with sign -1. 
#' At least one of 'lower' and 'higher' must be specified.
#' @param chunk_size number of database entries to process per iteration to 
#' limit memory usage of search.
#' @return \code{\link{gessResult}} object, the result table contains the 
#' search results for each perturbagen in the reference database ranked by 
#' their signature similarity to the query.
#' @seealso \code{\link{qSig}}, \code{\link{gessResult}}, \code{\link{gess}}
#' @references 
#' Lamb, J., Crawford, E. D., Peck, D., Modell, J. W., Blat, I. C., 
#' Wrobel, M. J., Golub, T. R. (2006). The Connectivity Map: 
#' using gene-expression signatures to connect small molecules, genes, and 
#' disease. Science, 313 (5795), 1929-1935. 
#' URL: https://doi.org/10.1126/science.1132939
#' 
#' Sandmann, T., Kummerfeld, S. K., Gentleman, R., & Bourgon, R. 
#' (2014). gCMAP: user-friendly connectivity mapping with R. Bioinformatics , 
#' 30 (1), 127-128. URL: https://doi.org/10.1093/bioinformatics/btt592
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