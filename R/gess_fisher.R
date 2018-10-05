#' Fisher method for GESS
#' 
#' @title gess_fisher
#' @param qSig `qSig` object, The 'gess_method' slot of 'qSig' should be 'Fisher'
#' @param higher The 'higher' threshold. If not 'NULL', genes with a score larger than 'higher' will be included in the gene set with sign +1. 
#' At least one of 'lower' and 'higher' must be specified.
#' @param lower The lower threshold. If not 'NULL', genes with a score smaller than 'lower' will be included in the gene set with sign -1. 
#' At least one of 'lower' and 'higher' must be specified.
#' @param chunk_size size of chunk per processing
#' @return gessResult object represents a list of drugs in reference database ranked by their similarity to query signature
#' @export
#' 
gess_fisher <- function(qSig, higher, lower, chunk_size=5000){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  #stopifnot(validObject(qSig))
  if(qSig@gess_method != "Fisher"){
    stop("The 'gess_method' slot of 'qSig' should be 'Fisher' if using 'gess_fisher' function")
  }
  query <- induceCMAPCollection(qSig@qsig, higher=higher, lower=lower)
  se <- qSig@refdb
  dmat <- assay(se)
  ceil <- ceiling(ncol(dmat)/chunk_size)
  resultDF <- data.frame()
  for(i in 1:ceil){
    dmat_sub <- dmat[,(chunk_size*(i-1)+1):min(chunk_size*i, ncol(dmat))]
    mat <- as(dmat_sub, "matrix")
    cmap <- induceCMAPCollection(mat, higher=higher, lower=lower)
    universe <- featureNames(cmap)
    c <- fisher_score(query=query, sets=cmap, universe = universe)
    resultDF <- rbind(resultDF, data.frame(c))
  }
  resultDF <- resultDF[order(resultDF$padj), ]
  row.names(resultDF) <- NULL
  x <- gessResult(result = as_tibble(resultDF),
                  qsig = qSig@qsig,
                  gess_method = qSig@gess_method,
                  refdb = qSig@refdb)
  return(x)
}