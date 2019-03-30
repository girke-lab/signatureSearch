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

##' @importFrom DOSE geneID
##' @method geneID feaResult
##' @export
geneID.feaResult <- function(x) as.character(x@result$geneID)

##' @importFrom DOSE geneInCategory
##' @importFrom stats setNames
##' @method geneInCategory feaResult
##' @export
geneInCategory.feaResult <- function(x)
  setNames(strsplit(geneID(x), "/", fixed=TRUE), x@result$ID)



