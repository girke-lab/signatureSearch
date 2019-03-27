#' @title Fisher method for GESS
#' @description 
#' It uses query signature to search against the reference database in the 
#' \code{qSig} by Fisher's exact test
#' @details 
#' The Fisher’s exact test can also be used as similarity search algorithm if 
#' both the query and the database are composed of DEG sets.
#' @param qSig `qSig` object, 
#' The 'gess_method' slot of 'qSig' should be 'Fisher'
#' @param higher The 'higher' threshold. If not 'NULL', genes with a score 
#' larger than 'higher' will be included in the gene set with sign +1. 
#' At least one of 'lower' and 'higher' must be specified.
#' @param lower The lower threshold. If not 'NULL', genes with a score smaller 
#' than 'lower' will be included in the gene set with sign -1. 
#' At least one of 'lower' and 'higher' must be specified.
#' @param chunk_size size of chunk per processing
#' @return gessResult object, containing drugs in the reference database
#' ranked by their similarity to the query signature
#' @seealso \code{\link{qSig}}, \code{\link{gessResult}}, \code{\link{gess}}
#' @examples 
#' db_dir <- system.file("extdata", "sample_db", package = "signatureSearch")
#' sample_db <- loadHDF5SummarizedExperiment(db_dir)
#' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
#' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
#' qsig_fisher <- qSig(qsig=query_mat, gess_method="Fisher", refdb=sample_db)
#' fisher <- gess_fisher(qSig=qsig_fisher, higher=1, lower=-1)
#' result(fisher)
#' @export
#' 
gess_fisher <- function(qSig, higher, lower, chunk_size=5000){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  #stopifnot(validObject(qSig))
  if(qSig@gess_method != "Fisher"){
    stop("The 'gess_method' slot of 'qSig' should be 'Fisher' 
         if using 'gess_fisher' function")
  }
  query <- induceCMAPCollection(qSig@qsig, higher=higher, lower=lower)
  se <- qSig@refdb
  dmat <- assay(se)
  ceil <- ceiling(ncol(dmat)/chunk_size)
  resultDF <- data.frame()
  for(i in seq_len(ceil)){
    dmat_sub <- dmat[,(chunk_size*(i-1)+1):min(chunk_size*i, ncol(dmat))]
    mat <- as(dmat_sub, "matrix")
    cmap <- induceCMAPCollection(mat, higher=higher, lower=lower)
    universe <- featureNames(cmap)
    c <- fisher_score(query=query, sets=cmap, universe = universe)
    resultDF <- rbind(resultDF, cmapTable(c))
  }
  resultDF <- resultDF[order(resultDF$padj), ]
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