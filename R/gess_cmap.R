#' @title CMAP Search Method
#' 
#' @description 
#' It uses query signature to search against the reference database defined
#' in the \code{\link{qSig}} object by CMAP method, which is an implementation 
#' of the CMap method from Lamb et al, 2006. 
#' 
#' @details 
#' Lamb et al. (2006) introduced the gene expression-based search method known 
#' as Connectivity Map (CMap) where a GES database is searched with a query GES 
#' for similar entries (Lamb et al. 2006). Specifically, the GESS method from 
#' Lamb et al. (2006), here termed as CMAP, uses as query the two label sets of 
#' the most up- and down-regulated genes from a genome-wide expression 
#' experiment, while the reference database is composed of rank transformed 
#' expression profiles (e.g. ranks of LFC or z-scores). The actual GESS 
#' algorithm is based on a vectorized rank difference calculation. The 
#' resulting Connectivity Score expresses to what degree the query up/down 
#' gene sets are enriched on the top and bottom of the database entries, 
#' respectively. The search results are a list of perturbagens such as drugs 
#' that induce similar or opposing GESs as the query. Similar GESs suggest 
#' similar physiological effects of the corresponding perturbagens. These GES 
#' associations can be useful to uncover novel MOAs of drugs or treatments for 
#' diseases. Although several variants of the CMAP algorithm are available in 
#' other software packages including Bioconductor, the implementation provided 
#' by this package follows the original description of the authors as closely 
#' as possible. This allows to reproduce in our tests the search results from 
#' the corresponding CMAP2 web service of the Broad Institute.
#' 
#' The CMAP algorithm takes about 1 minute on a single CPU core for querying 
#' with a single signature against 10,000 signatures in the reference database.
#' 
#' @section Column description:
#' Description of the score columns in the result table specific for CMAP 
#' method:
#' \itemize{
#'     \item raw_score: bi-directional enrichment score (Kolmogorov-Smirnov 
#'     statistic) of up and down set in the query siganture
#'     \item scaled_score: raw_score was scaled to valules from 1 to -1 by 
#'     dividing the positive scores with the maxmum positive score, and negative
#'     scores with the absolute value of minimum negative score.
#' }
#' Description of the other columns are available at the 'result' slot of the
#' \code{\link{gessResult}} object.
#' 
#' @param qSig \code{\link{qSig}} object defining the query signature, the GESS
#' method (should be 'CMAP') and the path to the reference database.
#' @param chunk_size size of chunk per processing
#' @return \code{\link{gessResult}} object, the result table contains the 
#' search results for each perturbagen in the reference database ranked by 
#' their signature similarity to the query.
#' @import methods
#' @seealso \code{\link{qSig}}, \code{\link{gessResult}}, \code{\link{gess}}
#' @references For detailed description of the CMap method, please refer to: 
#' Lamb, J., Crawford, E. D., Peck, D., Modell, J. W., Blat, I. C., 
#' Wrobel, M. J., … Golub, T. R. (2006). The Connectivity Map: 
#' using gene-expression signatures to connect small molecules, genes, and 
#' disease. Science, 313(5795), 1929–1935. 
#' \url{https://doi.org/10.1126/science.1132939}
#' @examples 
#' db_path <- system.file("extdata", "sample_db.h5", 
#'                        package = "signatureSearch")
#' library(signatureSearchData)
#' sample_db <- readHDF5chunk(db_path, colindex=1:100)
#' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
#' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
#' query = as.numeric(query_mat); names(query) = rownames(query_mat)
#' upset <- head(names(query[order(-query)]), 150)
#' downset <- tail(names(query[order(-query)]), 150)
#' qsig_cmap <- qSig(query = list(upset=upset, downset=downset), 
#'                   gess_method = "CMAP", refdb = db_path)
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
    db_path <- determine_refdb(qSig@refdb)
    qsig_up <- qSig@query[[1]]
    qsig_dn <- qSig@query[[2]]
    res <- cmapEnrich(db_path, upset=qsig_up, downset=qsig_dn, 
                      chunk_size=chunk_size)
    res <- sep_pcf(res)
    # add target column
    target <- suppressMessages(get_targets(res$pert))
    res <- left_join(res, target, by=c("pert"="drug_name"))
    
    x <- gessResult(result = as_tibble(res),
                    query = qSig@query,
                    gess_method = qSig@gess_method,
                    refdb = qSig@refdb)
    return(x)
}

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

#' @importFrom rhdf5 h5ls
cmapEnrich <- function(db_path, upset, downset, chunk_size=5000) {
  ## Read in matrix in h5 file by chunks
  mat_dim <- getH5dim(db_path)
  mat_nrow <- mat_dim[1]
  mat_ncol <- mat_dim[2]
  ceil <- ceiling(mat_ncol/chunk_size)
  ## get ranks of up and down genes in DB
  rankLup=NULL
  rankLdown=NULL
  for(i in seq_len(ceil)){
    mat <- readHDF5mat(db_path,
                    colindex=(chunk_size*(i-1)+1):min(chunk_size*i, mat_ncol))
    rankLup1 <- lapply(colnames(mat), function(x) sort(rank(-1*mat[,x])[upset]))
    rankLdown1 <- lapply(colnames(mat), 
                         function(x) sort(rank(-1*mat[,x])[downset]))
    rankLup <- c(rankLup, rankLup1)
    rankLdown <- c(rankLdown, rankLdown1)
  }
  
  ## Compute raw and scaled connectivity scores
  raw.score <- vapply(seq_along(rankLup), function(x) 
                               .s(rankLup[[x]], rankLdown[[x]], n=mat_nrow),
                      FUN.VALUE=numeric(1))
  score <- .S(raw.score)
  
  ## Assemble results
  resultDF <- data.frame(set = h5read(db_path, "colnames", drop=TRUE), 
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
