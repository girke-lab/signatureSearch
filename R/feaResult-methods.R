#' @name show
#' @docType methods
#' @rdname show-methods
#' @aliases show,feaResult-method
#' @importFrom utils str
#' @examples
#' fr <- feaResult(result=dplyr::tibble(id=letters[seq_len(10)], 
#'                                      val=seq_len(10)),
#'                 organism="human", ontology="MF", drugs=c("d1", "d2"), 
#'                 targets=c("t1","t2"))
#' fr 

setMethod("show", c(object="feaResult"),
      function (object){
          cat("#\n# Functional Enrichment Analysis \n#\n")
          cat("#...@organism", "\t", object@organism, "\n")
          cat("#...@ontology", "\t", object@ontology, "\n")
          cat("#...@drugs", "\t")
          str(object@drugs)
          cat("#...@targets", "\t")
          str(object@targets)
          # cat("#...@universe", "\t")
          # str(object@universe)
          cat(paste0("#...", nrow(result(object))), "enriched terms found\n")
          print(result(object))
          })

#' @name result
#' @docType methods
#' @rdname result-methods
#' @method result feaResult
#' @aliases result,feaResult-method
#' @examples
#' fr <- feaResult(result=dplyr::tibble(id=letters[seq_len(10)], 
#'                                      val=seq_len(10)),
#'                 organism="human", ontology="MF", drugs=c("d1", "d2"), 
#'                 targets=c("t1","t2"))
#' result(fr)

setMethod("result", c(x="feaResult"),
          function(x) x@result)

#' @description The \code{drugs} generic extracts or assign the drug names/ids 
#' stored in the \code{drugs} slot of an feaResult object.
#' @name drugs
#' @docType methods
#' @rdname drugs-methods
#' @method drugs feaResult
#' @aliases drugs,feaResult-method
#' @param x feaResult object
#' @return character vector
#' @examples 
#' fr <- feaResult(result=dplyr::tibble(id=letters[seq_len(10)], 
#'                                      val=seq_len(10)),
#'                 organism="human", ontology="MF", drugs=c("d1", "d2"), 
#'                 targets=c("t1","t2"))
#' drugs(fr)

setMethod("drugs", c(x="feaResult"),
          function(x) x@drugs)

#' @rdname drugs-methods
#' @method drugs feaResult
#' @aliases drugs,feaResult,ANY-method
#' @param value A character vector of drug names
#' @return An feaResult object with new assigned drugs slot
#' @examples 
#' drugs(fr) <- c("d3", "d4")
setMethod("drugs<-", "feaResult", function(x, value) {
    x@drugs <- value
    x
})

