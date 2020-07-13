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
#' larger than or equal to 'higher' will be included in the gene set with sign +1. 
#' At least one of 'lower' and 'higher' must be specified. 
#' 
#' \code{higher} argument need to be set as \code{1} if the \code{refdb} in \code{qSig} 
#' is path to the HDF5 file that were converted from the gmt file.
#' 
#' @param lower The lower threshold. If not 'NULL', genes with a score smaller 
#' than or equal 'lower' will be included in the gene set with sign -1. 
#' At least one of 'lower' and 'higher' must be specified.
#' 
#' \code{lower} argument need to be set as \code{NULL} if the \code{refdb} in \code{qSig} 
#' is path to the HDF5 file that were converted from the gmt file.
#' 
#' @param chunk_size number of database entries to process per iteration to 
#' limit memory usage of search.
#' @param ref_trts character vector. If users want to search against a subset 
#' of the reference database, they could set ref_trts as a character vector 
#' representing column names (treatments) of the subsetted refdb. 
#' @param workers integer(1) number of workers for searching the reference
#' database parallelly, default is 1.
#' 
#' @return \code{\link{gessResult}} object, the result table contains the 
#' search results for each perturbagen in the reference database ranked by 
#' their signature similarity to the query.
#' @importFrom HDF5Array HDF5Array
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
#' # library(SummarizedExperiment); library(HDF5Array)
#' # sample_db <- SummarizedExperiment(HDF5Array(db_path, name="assay"))
#' # rownames(sample_db) <- HDF5Array(db_path, name="rownames")
#' # colnames(sample_db) <- HDF5Array(db_path, name="colnames")
#' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
#' # query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
#' # qsig_gcmap <- qSig(query=query_mat, gess_method="gCMAP", refdb=db_path)
#' # gcmap <- gess_gcmap(qsig_gcmap, higher=1, lower=-1)
#' # result(gcmap)
#' @export
#' 
gess_gcmap <- function(qSig, higher=NULL, lower=NULL, chunk_size=5000, 
                       ref_trts=NULL, workers=1){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  #stopifnot(validObject(qSig))
  if(gm(qSig) != "gCMAP"){
    stop(paste("The 'gess_method' slot of 'qSig' should be 'gCMAP'",
               "if using 'gess_gcmap' function"))
  }
  if(!is.null(lower) && !is.null(higher) && higher==lower) {
    stop("Please specify two different cutoffs as 'higher' and 'lower'.")
  }
  query <- qr(qSig)
  db_path <- determine_refdb(refdb(qSig))
  
  ## calculate cs_raw of query to blocks (e.g., 5000 columns) of full refdb
  full_mat <- HDF5Array(db_path, "assay")
  rownames(full_mat) <- as.character(HDF5Array(db_path, "rownames"))
  colnames(full_mat) <- as.character(HDF5Array(db_path, "colnames"))
  
  if(! is.null(ref_trts)){
    trts_valid <- trts_check(ref_trts, colnames(full_mat))
    full_mat <- full_mat[, trts_valid]
  }
  
  full_dim <- dim(full_mat)
  full_grid <- colGrid(full_mat, ncol=min(chunk_size, ncol(full_mat)))
  ### The blocks in 'full_grid' are made of full columns 
  nblock <- length(full_grid) 
  resultDF <- bplapply(seq_len(nblock), function(b){
    ref_block <- read_block(full_mat, full_grid[[b]])
    cmap <- induceSMat(ref_block, higher=higher, lower=lower)
    c <- connectivity_score_raw(experiment=as.matrix(query), query=cmap)
    return(data.frame(c))}, BPPARAM = MulticoreParam(workers=workers))
  resultDF <- do.call(rbind, resultDF)

  # else {
  #   gene_sets <- suppressWarnings(read_gmt(db_path))
  #   # transform gene sets to sparseMatrix
  #   gsc <- GeneSetCollection(mapply(function(geneIds, setId) {
  #     GeneSet(geneIds, geneIdType=EntrezIdentifier(),
  #             setName=setId)
  #   }, gene_sets, names(gene_sets)))
  #   cmapData <- incidence(gsc)
  #   cmapData <- Matrix::t(cmapData)
  #   #c <- connectivity_score_raw(experiment=as.matrix(query), query=cmapData)
  #   # Search the refdb by using multiple cores/workers
  #   ceil <- ceiling(ncol(cmapData)/chunk_size)
  #   res_list <- bplapply(seq_len(ceil), function(i){
  #     dgcMat <- cmapData[,(chunk_size*(i-1)+1):min(chunk_size*i, ncol(cmapData))]
  #     c <- connectivity_score_raw(experiment=as.matrix(query), query=dgcMat)
  #     return(data.frame(c))
  #   }, BPPARAM = MulticoreParam(workers=workers))
  #   resultDF <- do.call(rbind, res_list)
  # }
  
  # mat_dim <- getH5dim(db_path)
  # mat_ncol <- mat_dim[2]
  # ceil <- ceiling(mat_ncol/chunk_size)
  # cs_raw <- function(i){
  #   mat <- readHDF5mat(db_path,
  #                   colindex=(chunk_size*(i-1)+1):min(chunk_size*i, mat_ncol))
  #   cmap <- gCMAP::induceCMAPCollection(mat, higher=higher, lower=lower)
  #   c <- connectivity_score_raw(experiment=as.matrix(query), query=cmap)
  #   return(data.frame(c))
  # }
  # resultDF <- do.call(rbind, lapply(seq_len(ceil), cs_raw))
  
  ## Apply scaling of scores to full data set
  resultDF[,"effect"] <-.connnectivity_scale(resultDF$effect)
  resultDF <- resultDF[order(abs(resultDF$effect), decreasing=TRUE), ]
  row.names(resultDF) <- NULL
  if(grepl("__.*__", resultDF$set[1])){
    resultDF <- sep_pcf(resultDF)
    # add target column
    target <- suppressMessages(get_targets(resultDF$pert))
    resultDF <- left_join(resultDF, target, by=c("pert"="drug_name"))
  }

  x <- gessResult(result = as_tibble(resultDF),
                  query = qr(qSig),
                  gess_method = gm(qSig),
                  refdb = refdb(qSig))
  return(x)
}

.connnectivity_scale <- function(scores) {
  p <- max(scores)
  q <- min(scores)
  ifelse(scores == 0, 0, ifelse( scores > 0, scores / p, -scores / q ))
}

#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix t
induceSMat <- function(mat, lower=NULL, higher=NULL){
  if( ! is.null(lower) && ! is.null(higher) && higher == lower) {
    stop("Please specify two different cutoffs")
  }
  gss <- lapply( seq_len(ncol(mat)),
                 function( n ) {
                   if (! is.null( lower )) {
                     down <- as.vector(which( mat[,n] <= lower ))
                   } else {
                     down <- NULL
                   }                            
                   if (! is.null( higher )) {
                     up <- as.vector(which( mat[,n] >= higher ))
                   } else {
                     up <- NULL
                   }                            
                   list( j = c(down, up),
                         x = c(rep(-1, length(down)), rep(1, length(up)))
                   )})
  i <- unlist(
    sapply( seq( length( gss ) ), function( m ) {
      rep( m, length( gss[[ m ]]$j ) )
    }))
  j <- unlist(sapply(gss ,function( m ) {m$j }))
  x <- unlist(sapply(gss ,function( m ) {m$x }))            
  cmap <- Matrix::t(sparseMatrix(i=as.integer(i),
                          j=as.integer(j),
                          x=as.integer(x),
                          dims=list(ncol(mat), nrow(mat)),
                          dimnames = list(colnames(mat), rownames(mat))))
  cmap
}

connectivity_score_raw <- function(experiment, query) {
  data.matrix <- experiment
  ## subset objects to shared genes
  matched.features <- fmatch(rownames(experiment), rownames(query))
  matched.sets <- query[na.omit(matched.features),]
  
  ## extract scores for each gene set
  sets.up <- lapply(seq(ncol(matched.sets)),
                    function(x) which(matched.sets[ ,x ] == 1))
  
  sets.down <- lapply(seq(ncol(matched.sets)),
                      function(x) which(matched.sets[ ,x] == -1))
  
  ## transform experiment to (reverse) ranks
  rank.matrix <- apply(data.matrix, 2, function(x) {length(x)-rank(x)+1})
  
  ## calculate connectivity score
  raw.score <- apply(rank.matrix, 2, function(r) {
    vapply(seq_along(sets.up), function(n) {
      .s(r[sets.up[[n]]], r[sets.down[[n]]], length(r))
    }, FUN.VALUE = numeric(1))
  })
  raw.score <- matrix(raw.score, ncol=ncol(experiment))
  
  ## scale score across all tested query sets
  ## ThG: modification on next two lines:
  ## score <- apply(raw.score, 2, .S)
  score <- raw.score
  score <- matrix(score, ncol=ncol(experiment))
  ## store results
  results <- data.frame(set = colnames(query), 
                        trend = ifelse(score[,1] >=0, "up", "down"),
                        effect = score[,1],
                        nSet = colSums(as.matrix(abs(query))),
                        nFound = colSums(as.matrix(abs(matched.sets))),
                        signed=TRUE)
}