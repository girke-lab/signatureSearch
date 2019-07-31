#' @title Fisher Search Method
#' @description 
#' It uses query signature to search against the reference database defined
#' in the \code{\link{qSig}} object by Fisher's exact test
#' @details 
#' Fisher’s exact test (Graham J. G. Upton 1992) can also be used to search a 
#' GES databases for entries that are similar to a query GES. In this case both 
#' the query and the database are composed of gene label sets, such as DEG sets.
#' 
#' @section Column description:
#' Description of the score columns in the result table specific for Fisher 
#' method:
#' \itemize{
#'     \item pval: p-value of the Fisher's exact test.
#'     \item padj: p-value adjusted for multiple hypothesis testing using
#'     R's p.adjust function with the Benjamini & Hochberg (BH) method. 
#'     \item effect: z-score based on the standard normal distribution.
#'     \item LOR: Log Odds Ratio.
#'     \item nSet: number of genes of the drug signature in the reference 
#'     database (gene sets) after setting the higher and lower cutoff.
#'     \item nFound: number of genes of the drug signature in the reference 
#'     database (gene sets) also found in the query signature
#'     (whole genome profile).
#'     \item signed: whether gene sets in the reference database have signs, 
#'     representing up and down regulated genes when computing scores.  
#' }
#' Description of the other columns are available at the 'result' slot of the
#' \code{\link{gessResult}} object.
#' 
#' @param qSig \code{\link{qSig}} object defining the query signature, the GESS
#' method (should be 'Fisher') and the path to the reference database.
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
#' @importMethodsFrom GSEABase GeneSet
#' @seealso \code{\link{qSig}}, \code{\link{gessResult}}, \code{\link{gess}}
#' @references 
#' Graham J. G. Upton. 1992. “Fisher’s Exact Test.” J. R. Stat. Soc. Ser. A 
#' Stat. Soc. 155 (3). [Wiley, Royal Statistical Society]: 395–402. 
#' \url{http://www.jstor.org/stable/2982890}
#' @examples 
#' db_path <- system.file("extdata", "sample_db.h5", 
#'                        package = "signatureSearch")
#' library(signatureSearchData)
#' sample_db <- readHDF5chunk(db_path, colindex=1:100)
#' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
#' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
#' qsig_fisher <- qSig(query=query_mat, gess_method="Fisher", refdb=db_path)
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
  if(is.list(qSig@query)){
    query <- GeneSet(unique(unlist(qSig@query)))
  } else {
    query <- induceCMAPCollection(qSig@query, higher=higher, lower=lower)
  }
  db_path <- determine_refdb(qSig@refdb)
  mat_dim <- getH5dim(db_path)
  mat_ncol <- mat_dim[2]
  ceil <- ceiling(mat_ncol/chunk_size)
  resultDF <- data.frame()
  for(i in seq_len(ceil)){
    mat <- readHDF5mat(db_path,
                    colindex=(chunk_size*(i-1)+1):min(chunk_size*i, mat_ncol))
    cmap <- induceCMAPCollection(mat, higher=higher, lower=lower)
    universe <- featureNames(cmap)
    c <- fisher_score(query=query, sets=cmap, universe = universe)
    resultDF <- rbind(resultDF, cmapTable(c))
  }
  resultDF <- resultDF[order(resultDF$padj), ]
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