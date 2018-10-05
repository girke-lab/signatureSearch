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

##' @method geneID feaResult
##' @export
geneID.feaResult <- function(x) as.character(x@result$geneID)


##' @method geneInCategory feaResult
##' @export
##' @importFrom stats setNames
geneInCategory.feaResult <- function(x)
    setNames(strsplit(geneID(x), "/", fixed=TRUE), rownames(x@result))

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

##' get 'result' slot of gessResult object
##' @name result
##' @docType methods
##' @rdname result-methods
##' @method result gessResult
##' @param x gessResult object
##' @return tibble
##' @aliases result,gessResult-method
setMethod("result", signature(x="gessResult"),
          function(x) x@result)

## Constructor for "gessResult"
gessResult <- function(result, qsig, gess_method, refdb)
  new("gessResult", result=result, qsig=qsig, gess_method=gess_method, refdb=refdb)

##' show qSig, gessResult, feaResult objects
##' 
##' @name show
##' @docType methods
##' @rdname show-methods
##' @title show method
##' @param object object used for show
##' @return message
##' @aliases show,gessResult-method
##' @usage show(object)
setMethod("show", signature(object="gessResult"),
          function (object) {
            cat("#\n# gessResult object \n#\n")
            cat("@result \n")
            print(object@result)
            if(is(object@qsig, "list")){
              cat("@qsig", "\t", "up gene set", paste0("(", length(object@qsig[[1]]), "):"), "\t", object@qsig[[1]][1:10], "... \n")
              cat("     ", "\t", "down gene set", paste0("(", length(object@qsig[[2]]), "):"), "\t", object@qsig[[2]][1:10], "... \n")
            }
            if(is(object@qsig, "matrix")){
              cat("@qsig\n")
              mat=object@qsig
              print(head(mat,10))
              cat("# ... with", nrow(mat)-10, "more rows\n")
            }
            cat("\n@gess_method", "\t", object@gess_method, "\n")
            cat("\n@refdb", "\t")
            print(object@refdb)
          })



