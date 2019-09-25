##' Return the first part of the result table in the \code{\link{gessResult}}, 
##' and \code{\link{feaResult}} objects
##' @title Return the First Part of an Object 
##' @importFrom utils head
##' @name head
##' @docType methods
##' @rdname head-methods
##' @param x an object
##' @param n a single integer. If positive or zero, size for the resulting 
##' object is the number of rows for a data frame. If negative, all but the n 
##' last number of rows of x.
##' @param ... arguments to be passed to or from other methods
##' @return data.frame
##' @aliases head,gessResult-method
##' @method head gessResult
##' @examples 
##' gr <- gessResult(result=dplyr::tibble(pert=letters[seq_len(10)], 
##'                                       val=seq_len(10)), 
##'                  query=list(up=c("g1","g2"), down=c("g3","g4")),
##'                  gess_method="LINCS", refdb="path/to/lincs/db")
##' head(gr)
##' @export
setMethod("head", "gessResult",
          function(x, n=6L, ...) utils::head(as.data.frame(x@result), n, ...))

##' @name head
##' @docType methods
##' @rdname head-methods
##' @aliases head,feaResult-method
##' @method head feaResult
##' @examples 
##' fr <- feaResult(result=dplyr::tibble(id=letters[seq_len(10)], 
##'                                     val=seq_len(10)),
##'                 organism="human", ontology="MF", drugs=c("d1", "d2"), 
##'                 targets=c("t1","t2"))
##' head(fr)
##' @export
setMethod("head", "feaResult",
          function(x, n=6L, ...) utils::head(as.data.frame(x@result), n, ...))

##' Return the last part of the result table in the \code{\link{gessResult}}, 
##' and \code{\link{feaResult}} objects
##' @title Return the Last Part of an Object 
##' @importFrom utils tail
##' @name tail
##' @docType methods
##' @rdname tail-methods
##' @param x an object
##' @param n a single integer. If positive or zero, size for the resulting 
##' object is the number of rows for a data frame. If negative, all but the n 
##' first number of rows of x.
##' @param ... arguments to be passed to or from other methods
##' @return data.frame
##' @aliases tail,gessResult-method
##' @method tail gessResult
##' @examples 
##' gr <- gessResult(result=dplyr::tibble(pert=letters[seq_len(10)], 
##'                                       val=seq_len(10)), 
##'                  query=list(up=c("g1","g2"), down=c("g3","g4")),
##'                  gess_method="LINCS", refdb="path/to/lincs/db")
##' tail(gr)
##' @export
setMethod("tail", "gessResult",
          function(x, n=6L, ...) utils::tail(as.data.frame(x@result), n, ...))

##' @name tail
##' @docType methods
##' @rdname tail-methods
##' @aliases tail,feaResult-method
##' @method tail feaResult
##' @examples 
##' fr <- feaResult(result=dplyr::tibble(id=letters[seq_len(10)], 
##'                                      val=seq_len(10)),
##'                 organism="human", ontology="MF", drugs=c("d1", "d2"), 
##'                 targets=c("t1","t2"))
##' tail(fr)
##' @export
setMethod("tail", "feaResult",
          function(x, n=6L, ...) utils::tail(as.data.frame(x@result), n, ...))

##' Retrieve dimension of the result table in the \code{\link{gessResult}}, 
##' and \code{\link{feaResult}} objects
##' @title Dimensions of an Object
##' @name dim
##' @docType methods
##' @rdname dim-methods
##' @param x an R object
##' @return dim attribute of the result table
##' @aliases dim,gessResult-method
##' @method dim gessResult
##' @examples 
##' gr <- gessResult(result=dplyr::tibble(pert=letters[seq_len(10)], 
##'                                       val=seq_len(10)), 
##'                  query=list(up=c("g1","g2"), down=c("g3","g4")),
##'                  gess_method="LINCS", refdb="path/to/lincs/db")
##' dim(gr)
##' @export
setMethod("dim", "gessResult",
          function(x) dim(x@result))

##' @name dim
##' @docType methods
##' @rdname dim-methods
##' @aliases dim,feaResult-method
##' @method dim feaResult
##' @examples 
##' fr <- feaResult(result=dplyr::tibble(id=letters[seq_len(10)], 
##'                                      val=seq_len(10)),
##'                 organism="human", ontology="MF", drugs=c("d1", "d2"), 
##'                 targets=c("t1","t2"))
##' dim(fr)
##' @export
setMethod("dim", "feaResult",
          function(x) dim(x@result))

# setMethod("[", "feaResult",
#           function(x, i, j) x@result[i,j])
`[.feaResult` <- function(x, i, j) {
    x@result[i,j]
}

# setMethod("$", "feaResult",
#           function(x, name) x@result[, name])
`$.feaResult` <-  function(x, name) {
    x@result[, name]
}

qr <- function(x) x@query
gm <- function(x) x@gess_method
refdb <- function(x) x@refdb

tg <- function(x) x@targets
`tg<-` <- function(x, value){
    x@targets <- value
    return(x)}
og <- function(x) x@organism
`og<-` <- function(x, value){
    x@organism <- value
    return(x)}
ont <- function(x) x@ontology
`ont<-` <- function(x, value){
    x@ontology <- value
    return(x)}
`rst<-` <- function(x, value){
    x@result <- value
    return(x)}


