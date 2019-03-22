#' @importFrom data.table frank
rankMatrix <- function(x, decreasing=TRUE) {
  if(!is(x, "matrix") | !is.numeric(x)) stop("x needs to be numeric matrix!")
  if(is.null(rownames(x)) | is.null(colnames(x))) 
    stop("matrix x lacks colnames and/or rownames!")
  if(decreasing==TRUE) { mysign <- -1 } else { mysign <- 1 }
  rankma <- vapply(seq(ncol(x)), function(z) data.table::frank(mysign * x[,z]))
  rownames(rankma) <- rownames(x); colnames(rankma) <- colnames(x)
  return(rankma)
}

cmapEnrich <- function(se, upset, downset, chunk_size=5000) {
  ## Validity checks of inputs
  # if(!is(rankMA, "matrix") | !is.numeric(rankMA)) 
  #     stop("rankMA needs to be numeric matrix!")
  # if(is.null(rownames(rankMA)) | is.null(colnames(rankMA))) 
  #     stop("matrix rankMA lacks colnames and/or rownames!")
  dmat <- assay(se)
  ## Obtain ranks of query genes in DB
  ceil <- ceiling(ncol(dmat)/chunk_size)
  rankLup=NULL
  rankLdown=NULL
  for(i in seq_len(ceil)){
    dmat_sub <- dmat[,(chunk_size*(i-1)+1):min(chunk_size*i, ncol(dmat))]
    mat <- as(dmat_sub, "matrix")
    rankLup1 <- lapply(colnames(mat), function(x) sort(rank(-1*mat[,x])[upset]))
    rankLdown1 <- lapply(colnames(mat), 
                         function(x) sort(rank(-1*mat[,x])[downset]))
    rankLup <- c(rankLup, rankLup1)
    rankLdown <- c(rankLdown, rankLdown1)
  }
  
  ## Compute raw and scaled connectivity scores
  raw.score <- vapply(seq_along(rankLup), function(x) 
                               .s(rankLup[[x]], rankLdown[[x]], n=nrow(dmat)),
                      FUN.VALUE=numeric(1))
  score <- .S(raw.score)
  
  ## Assemble results
  resultDF <- data.frame(set = colnames(dmat), 
                         trend = ifelse(score >=0, "up", "down"),
                         raw_score = raw.score,
                         scaled_score = score,
                         N_upset = length(upset),
                         N_downset = length(downset), stringsAsFactors = FALSE)
  resultDF <- resultDF[order(abs(resultDF$scaled_score), decreasing=TRUE), ]
  rownames(resultDF) <- NULL
  return(resultDF)
}

## Fct to compute a and b
.ks <- function(V, n) {
  t <- length(V)
  a <- max(seq_len(t) / t - V / n)
  b <- max(V / n - (seq_len(t)-1)/t)
  ifelse(a > b, a, -b)
}

## Fct to compute ks_up and ks_down
.s <- function(V_up, V_down, n) {
  ks_up <- .ks(V_up, n)
  ks_down <- .ks(V_down, n)
  ifelse(sign(ks_up) == sign(ks_down), 0, ks_up - ks_down)
}

## Fct to scale scores
.S <- function(scores) {
  p <- max(scores)
  q <- min(scores)
  ifelse(scores == 0, 0, ifelse(scores > 0, scores / p, -scores / q))
}

#' @title CMAP method for GESS
#' 
#' @description 
#' It uses query signature to search against the reference database in the 
#' \code{qSig} by CMAP method, which is an implementation of the original 
#' CMap method from Lamb et al, 2006. 
#' 
#' @details 
#' Lamb et at., 2006 introduced the gene expression-based search method known as 
#' Connectivity Map (CMap) where a GES database is searched with a query GES for 
#' similar matches. It uses as query the most up- and down-regulated DEGs from 
#' a genome-wide expression experiment. The GES query is used to search a 
#' database of rank transformed GEPs and ranks the results by the degree of 
#' enrichment of the up- and down-regulated query genes on the top and bottom 
#' of the databases entries, respectively. The search results are a list of 
#' drugs that have similar or opposing GESs as the query. Similar GESs suggest 
#' similar physiological effects of the corresponding perturbagens. 
#' 
#' The CMAP method takes about 4 min on a single CPU core for querying with a 
#' single signature against ~100,000 signatures in the LINCS database.
#' 
#' @param qSig \code{qSig} object, The 'gess_method' slot should be 'CMAP'.
#' @param chunk_size size of chunk per processing
#' @return \code{gessResult} object, represents drugs in the reference database 
#' ranked by their similarity to the query signature
#' @importFrom utils download.file
#' @importFrom R.utils gunzip
#' @importFrom tools file_path_as_absolute
#' @import methods
#' @seealso \code{\link{qSig}}, \code{\link{gessResult}}, \code{\link{gess}}
#' @references For detailed description of the CMap method, please refer to
#' Lamb et al., 2006, \url{http://science.sciencemag.org/content/313/5795/1929}
#' @examples 
#' db_dir <- system.file("extdata", "sample_db", package = "signatureSearch")
#' sample_db <- loadHDF5SummarizedExperiment(db_dir)
#' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
#' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
#' query = as.numeric(query_mat); names(query) = rownames(query_mat)
#' upset <- head(names(query[order(-query)]), 150)
#' downset <- tail(names(query[order(-query)]), 150)
#' qsig_cmap <- qSig(qsig = list(upset=upset, downset=downset), 
#'                   gess_method = "CMAP", refdb = sample_db)
#' cmap <- gess_cmap(qSig=qsig_cmap, chunk_size=5000)
#' result(cmap)
#' @export
gess_cmap <- function(qSig, chunk_size=5000){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  # stopifnot(validObject(qSig))
  if(qSig@gess_method != "CMAP"){
    stop(paste("The 'gess_method' slot of 'qSig' should be 'CMAP'",
               "if using 'gess_cmap' function!"))
  }
  se <- qSig@refdb
  #dmat = assay(se)
  qsig_up <- qSig@qsig[[1]]
  qsig_dn <- qSig@qsig[[2]]
  res <- cmapEnrich(se, upset=qsig_up, downset=qsig_dn, chunk_size=chunk_size)
  new <- as.data.frame(t(vapply(seq_len(nrow(res)), function(i)
    unlist(strsplit(as.character(res$set[i]), "__")), FUN.VALUE=character(3))), 
    stringsAsFactors=FALSE)
  colnames(new) = c("pert", "cell", "type")
  res <- cbind(new, res[,-1])
  # add target column
  target <- suppressMessages(get_targets(res$pert))
  res <- left_join(res, target, by=c("pert"="drug_name"))
    
  x <- gessResult(result = as_tibble(res),
                  qsig = qSig@qsig,
                  gess_method = qSig@gess_method,
                  refdb = qSig@refdb)
  return(x)
}
