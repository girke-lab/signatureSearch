#' Rank transform signature database provided as matrix
#' @title rankMatrix
#' @param x GEPs
#' @param decreasing rank decreasingly
#' @importFrom data.table frank
#' @examples 
#' \dontrun{
#'  utils::download.file("http://biocluster.ucr.edu/~yduan004/CMAP_db/degList.rds", 
#'  "degList.rds", quiet = TRUE)
#'  logMA <- readRDS("degList.rds")$logFC
#   rankMA <- rankMatrix(x=logMA, decreasing=TRUE)
#' }
rankMatrix <- function(x, decreasing=TRUE) {
  if(class(x)!="matrix" | !is.numeric(x)) stop("x needs to be numeric matrix!")
  if(is.null(rownames(x)) | is.null(colnames(x))) stop("matrix x lacks colnames and/or rownames!")
  if(decreasing==TRUE) { mysign <- -1 } else { mysign <- 1 }
  rankma <- sapply(seq(ncol(x)), function(z) data.table::frank(mysign * x[,z]))
  rownames(rankma) <- rownames(x); colnames(rankma) <- colnames(x)
  return(rankma)
}

cmapEnrich <- function(se, upset, downset, chunk_size=5000) {
  ## Validity checks of inputs
  # if(class(rankMA)!="matrix" | !is.numeric(rankMA)) stop("rankMA needs to be numeric matrix!")
  # if(is.null(rownames(rankMA)) | is.null(colnames(rankMA))) stop("matrix rankMA lacks colnames and/or rownames!")
  
  dmat <- assay(se)
  ## Obtain ranks of query genes in DB
  ceil <- ceiling(ncol(dmat)/chunk_size)
  rankLup=NULL
  rankLdown=NULL
  for(i in 1:ceil){
    dmat_sub <- dmat[,(chunk_size*(i-1)+1):min(chunk_size*i, ncol(dmat))]
    mat <- as(dmat_sub, "matrix")
    rankLup1 <- sapply(colnames(mat), function(x) sort(rank(-1*mat[,x])[upset]), simplify=FALSE)
    rankLdown1 <- sapply(colnames(mat), function(x) sort(rank(-1*mat[,x])[downset]), simplify=FALSE)
    rankLup <- c(rankLup, rankLup1)
    rankLdown <- c(rankLdown, rankLdown1)
  }
  
  ## Fct to compute a and b
  .ks <- function(V, n) {
    t <- length(V)
    a <- max((1:t) / t - V / n)
    b <- max(V / n - ((1:t)-1)/t)
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
  
  ## Compute raw and scaled connectivity scores
  raw.score <- sapply(seq_along(rankLup), function(x) .s(rankLup[[x]], rankLdown[[x]], n=nrow(dmat)))
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

#' CMAP method for GESS
#' 
#' Motivation: gene expression signature search (GESS) using original CMap enrichment and new LINCS algorithms, as well as variants from gCMAP package
#' against very large expression databases containing tens to hundreds of thousands of profiles with very moderate memory requirements. The original 
#' CMap and new LINCS enrichment functions are custom implementations according to Lamb et al, 2006 and Subramanian et al, 2017. Most methods 
#' can be easily parallelized for multiple query signatures. The time performance of the gess_cmap/gess_lincs methods on a single CPU core is about 
#' 4 min for querying with a single signature ~100,000 CMAP/LINCS database signatures for gess_cmap method, while it is ~10 min for gess_lincs method.
#' 
#' @title gess_cmap
#' @param qSig `qSig` object, The 'gess_method' slot should be 'CMAP'. The GEPs will be internally rank transformed.
#' @param chunk_size size of chunk per processing
#' @return `gessResult` object, represents a list of drugs in reference database ranked by their similarity to query signature
#' @importFrom utils download.file
#' @importFrom R.utils gunzip
#' @importFrom tools file_path_as_absolute
#' @import methods
#' @export
#' 
gess_cmap <- function(qSig, chunk_size=5000){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  # stopifnot(validObject(qSig))
  if(qSig@gess_method != "CMAP"){
    stop("The 'gess_method' slot of 'qSig' should be 'CMAP' if using 'gess_cmap' function!")
  }
  se <- qSig@refdb
  #dmat = assay(se)
  qsig_up <- qSig@qsig[[1]]
  qsig_dn <- qSig@qsig[[2]]
  res <- cmapEnrich(se, upset=qsig_up, downset=qsig_dn, chunk_size=chunk_size)
  new <- as.data.frame(t(sapply(1:nrow(res), function(i)
    unlist(strsplit(as.character(res$set[i]), "__")))), stringsAsFactors=FALSE)
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
