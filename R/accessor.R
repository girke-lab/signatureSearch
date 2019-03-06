##' @method as.data.frame feaResult
##' @export
as.data.frame.feaResult <- function(x, ...) {
    as.data.frame(x@result, ...)
}

##' @method as.data.frame gessResult
##' @export
as.data.frame.gessResult <- function(x, ...) {
  as.data.frame(x@result, ...)
}

##' @method [ feaResult
##' @export
`[.feaResult` <- function(x, i, j) {
  x@result[i,j]
}

##' @method [ gessResult
##' @export
`[.gessResult` <- function(x, i, j) {
  x@result[i,j]
}

##' @method $ feaResult
##' @export
`$.feaResult` <-  function(x, name) {
    x@result[, name]
}

##' @method $ gessResult
##' @export
`$.gessResult` <-  function(x, name) {
  x@result[, name]
}

##' @method [[ feaResult
##' @export
`[[.feaResult` <- function(x, i) {
    gc <- geneInCategory(x)
    if (!i %in% names(gc))
        stop("input term not found...")
    gc[[i]]
}

##' @importFrom utils head
##' @method head feaResult
##' @export
head.feaResult <- function(x, n=6L, ...) {
    utils::head(x@result, n, ...)
}

##' @importFrom utils head
##' @method head gessResult
##' @export
head.gessResult <- function(x, n=6L, ...) {
  utils::head(x@result, n, ...)
}

##' @importFrom utils tail
##' @method tail feaResult
##' @export
tail.feaResult <- function(x, n=6L, ...) {
    utils::tail(x@result, n, ...)
}

##' @importFrom utils tail
##' @method tail gessResult
##' @export
tail.gessResult <- function(x, n=6L, ...) {
  utils::tail(x@result, n, ...)
}

##' @method dim feaResult
##' @export
dim.feaResult <- function(x) {
    dim(x@result)
}

##' @method dim gessResult
##' @export
dim.gessResult <- function(x) {
  dim(x@result)
}

##' geneID generic
##' @name geneID
##' @rdname geneID-method
##' @method geneID feaResult
##' @docType methods
##' @param x feaResult object
##' @return 'geneID' return the 'geneID' column of the FEA result which can be converted to data.frame via 'as.data.frame'
##' @importFrom DOSE geneID
##' @export
##' @examples
##' data(drugs, package="signatureSearch")
##' dup_hyperG_res <- tsea_dup_hyperG(drugs = drugs, universe = "Default", type = "GO", ont="MF")
##' head(geneID(dup_hyperG_res))
geneID.feaResult <- function(x) as.character(x@result$geneID)

##' geneInCategory generic
##' @name geneInCategory
##' @rdname geneInCategory-method
##' @method geneInCategory feaResult
##' @docType methods
##' @param x feaResult
##' @return 'geneInCategory' return a list of genes, by spliting the input gene vector to enriched functional categories
##' @importFrom stats setNames
##' @importFrom DOSE geneInCategory
##' @export
##' @examples
##' data(drugs, package="signatureSearch")
##' dup_hyperG_res <- tsea_dup_hyperG(drugs = drugs, universe = "Default", type = "GO", ont="MF")
##' head(geneInCategory(dup_hyperG_res))
geneInCategory.feaResult <- function(x)
  setNames(strsplit(geneID(x), "/", fixed=TRUE), x@result$ID)



