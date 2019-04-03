#' @title gCMAP method for GESS
#' @description 
#' It uses query signature to search against the reference database in the 
#' \code{qSig} by gCMAP method, which is adapted from the gCMAP package 
#' (Sandmann et al., 2014)
#' @details 
#' The \code{gCMAP} package provides access to related 
#' but not identical implementations of the original CMAP algorithm proposed by 
#' Lamb et al., 2006. It uses as query a rank transformed GEP and the reference 
#' database is composed of DEG sets. This is the opposite situation of the 
#' original \code{CMAP} method, where the query is a DEG set and the database 
#' contains rank transformed GEPs.
#' 
#' Description of the score columns in the gess_gcmap tibble result:
#' \itemize{
#'     \item effect: scaled bi-directional enrichment score, the same as 
#'     the 'scaled_score' in the \code{\link{gess_cmap}} result
#'     \item nSet: number of genes of the drug signature in the reference 
#'     database (gene sets) after setting the higher and lower cutoff.
#'     \item nFound: number of genes of the drug signature in the reference 
#'     database (gene sets) also found in the query signature 
#'     (whole genome profile).
#'     \item signed: whether gene sets in the reference database have signs, 
#'     representing up and down regulated genes when computing scores 
#' }
#' @param qSig `qSig` object, The 'gess_method' slot of 'qSig' should be 'gCMAP'
#' @param higher The 'higher' threshold. If not 'NULL', genes with a score 
#' larger than 'higher' will be included in the gene set with sign +1. 
#' At least one of 'lower' and 'higher' must be specified.
#' @param lower The lower threshold. If not 'NULL', genes with a score smaller 
#' than 'lower' will be included in the gene set with sign -1. 
#' At least one of 'lower' and 'higher' must be specified.
#' @param add_bs_score TRUE or FALSE. If true, bootstrap scores are added to 
#' measure the robustness of the result rankings.
#' @param chunk_size size of chunk per processing
#' @return gessResult object, containing drugs in the reference database
#' ranked by their similarity to the query signature
#' @importFrom R.utils gunzip
#' @import HDF5Array
#' @import SummarizedExperiment
#' @seealso \code{\link{qSig}}, \code{\link{gessResult}}, \code{\link{gess}}
#' @references 
#' \itemize{
#'   \item Sandmann et al., 2014,
#'         \url{https://academic.oup.com/bioinformatics/article/30/1/127/236809}
#'   \item gCMAP package,
#'         \url{https://bioconductor.org/packages/release/bioc/html/gCMAP.html}
#' }
#' @examples 
#' db_dir <- system.file("extdata", "sample_db", package = "signatureSearch")
#' sample_db <- loadHDF5SummarizedExperiment(db_dir)
#' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
#' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
#' qsig_gcmap <- qSig(qsig=query_mat, gess_method="gCMAP", refdb=sample_db,
#'                    refdb_name="sample")
#' gcmap <- gess_gcmap(qsig_gcmap, higher=1, lower=-1)
#' result(gcmap)
#' @export
#' 
gess_gcmap <- function(qSig, higher, lower, 
                       add_bs_score = FALSE, chunk_size=5000){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  #stopifnot(validObject(qSig))
  if(qSig@gess_method != "gCMAP"){
    stop(paste("The 'gess_method' slot of 'qSig' should be 'gCMAP'",
               "if using 'gess_gcmap' function"))
  }
  query <- qSig@qsig
  se <- qSig@refdb
  dmat <- assay(se)
  ceil <- ceiling(ncol(dmat)/chunk_size)
  resultDF <- data.frame()
  for(i in seq_len(ceil)){
    dmat_sub <- dmat[,(chunk_size*(i-1)+1):min(chunk_size*i, ncol(dmat))]
    mat <- as(dmat_sub, "matrix")
    cmap <- gCMAP::induceCMAPCollection(mat, higher=higher, lower=lower)
    c <- connectivity_score_raw(experiment=as.matrix(query), query=cmap)
    resultDF <- rbind(resultDF, data.frame(c))
  }
  ## Apply scaling of scores to full data set
  resultDF[,"effect"] <-.connnectivity_scale(resultDF$effect)
  resultDF <- resultDF[order(abs(resultDF$effect), decreasing=TRUE), ]
  row.names(resultDF) <- NULL
  
  new <- as.data.frame(t(vapply(seq_len(nrow(resultDF)), function(i)
    unlist(strsplit(as.character(resultDF$set[i]), "__")),
    FUN.VALUE = character(3))), stringsAsFactors=FALSE)
  colnames(new) = c("pert", "cell", "type")
  resultDF <- cbind(new, resultDF[,-1])
  # add target column
  target <- suppressMessages(get_targets(resultDF$pert))
  res <- left_join(resultDF, target, by=c("pert"="drug_name"))
  
  x <- gessResult(result = as_tibble(res),
                  qsig = qSig@qsig,
                  gess_method = qSig@gess_method,
                  refdb_name = qSig@refdb_name)
  return(x)
}

.connnectivity_scale <- function(scores) {
  p <- max(scores)
  q <- min(scores)
  ifelse(scores == 0, 0, ifelse( scores > 0, scores / p, -scores / q ))
}