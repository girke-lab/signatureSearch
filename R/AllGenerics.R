##' Generate \code{qSig} object used for GESS methods
##' 
##' @title qSig method
##' @rdname qSig-methods
##' @export
setGeneric("qSig", function(query, gess_method, refdb) 
  standardGeneric("qSig"))

##' Get GESS or FEA result tables
##' 
##' @title result method
##' @rdname result-methods
##' @export
setGeneric("result", function(x) standardGeneric("result"))

##' get_drugs generic
##' 
##' @title get_drugs method
##' @rdname get_drugs-methods
##' @export
setGeneric("get_drugs", function(x) standardGeneric("get_drugs"))

setGeneric("connectivity_score_raw",
           function( experiment, query, ... ){
             standardGeneric( "connectivity_score_raw" )
})
