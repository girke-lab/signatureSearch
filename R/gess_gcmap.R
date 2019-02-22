#' gCMAP method for GESS
#' 
#' @title gess_gcmap
#' @param qSig `qSig` object, The 'gess_method' slot of 'qSig' should be 'gCMAP'
#' @param higher The 'higher' threshold. If not 'NULL', genes with a score larger than 'higher' will be included in the gene set with sign +1. 
#' At least one of 'lower' and 'higher' must be specified.
#' @param lower The lower threshold. If not 'NULL', genes with a score smaller than 'lower' will be included in the gene set with sign -1. 
#' At least one of 'lower' and 'higher' must be specified.
#' @param add_bs_score TRUE or FALSE. If true, bootstrap scores are added to measure the robustness of the result rankings.
#' @param chunk_size size of chunk per processing
#' @return gessResult object represents a list of drugs in reference database ranked by their similarity to query signature
#' @importFrom utils download.file
#' @importFrom R.utils gunzip
#' @import HDF5Array
#' @import SummarizedExperiment
#' @export
#' 
gess_gcmap <- function(qSig, higher, lower, add_bs_score = FALSE, chunk_size=5000){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  #stopifnot(validObject(qSig))
  if(qSig@gess_method != "gCMAP"){
    stop("The 'gess_method' slot of 'qSig' should be 'gCMAP' if using 'gess_gcmap' function")
  }
  query <- qSig@qsig
  se <- qSig@refdb
  dmat <- assay(se)
  ceil <- ceiling(ncol(dmat)/chunk_size)
  resultDF <- data.frame()
  for(i in 1:ceil){
    dmat_sub <- dmat[,(chunk_size*(i-1)+1):min(chunk_size*i, ncol(dmat))]
    mat <- as(dmat_sub, "matrix")
    cmap <- induceCMAPCollection(mat, higher=higher, lower=lower)
    c <- connectivity_score(experiment=as.matrix(query), query=cmap)
    resultDF <- rbind(resultDF, data.frame(c))
  }
  ## Apply scaling of scores to full data set
  .connnectivity_scale <- function(scores) {
    p <- max(scores)
    q <- min(scores)
    ifelse(scores == 0, 0, ifelse( scores > 0, scores / p, -scores / q ))
  }
  resultDF[,"effect"] <-.connnectivity_scale(resultDF$effect)
  resultDF <- resultDF[order(abs(resultDF$effect), decreasing=TRUE), ]
  row.names(resultDF) <- NULL
  
  new <- as.data.frame(t(sapply(1:nrow(resultDF), function(i)
    unlist(strsplit(as.character(resultDF$set[i]), "__")))), stringsAsFactors=FALSE)
  colnames(new) = c("pert", "cell", "type")
  resultDF <- cbind(new, resultDF[,-1])
  # add target column
  target <- suppressMessages(get_targets(resultDF$pert))
  res <- left_join(resultDF, target, by=c("pert"="drug_name"))
  
  x <- gessResult(result = as_tibble(res),
                  qsig = qSig@qsig,
                  gess_method = qSig@gess_method,
                  refdb = qSig@refdb)
  return(x)
}