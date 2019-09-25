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
          cat(paste0("#...", nrow(object@result)), "enriched terms found\n")
          str(object@result)
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

# ##' dotplot for feaResult
# ##'
# ##' @rdname dotplot-methods
# ##' @aliases dotplot,feaResult,ANY-method
# ##' @param object an instance of feaResult
# ##' @param x variable for x axis
# ##' @param colorBy one of 'pvalue', 'p.adjust' and 'qvalue'
# ##' @param showCategory number of category
# ##' @param split separate result by 'category' variable
# ##' @param font.size font size
# ##' @param title plot title
# ##' @return plot
# ##' @importFrom enrichplot dotplot
# ##' @exportMethod dotplot
# setMethod("dotplot", signature(object="feaResult"),
#           function(object, x="geneRatio", colorBy="p.adjust", 
#           showCategory=10, split=NULL, font.size=12, title="") {
#               dotplot_internal(object, x, colorBy, showCategory, 
#                                split, font.size, title)
#            })

