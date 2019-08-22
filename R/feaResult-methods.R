##' @name show
##' @docType methods
##' @rdname show-methods
##' @aliases show,feaResult-method
##' @importFrom utils str
##' @examples
##' data(drugs)
##' # dup_hyperG_res <- tsea_dup_hyperG(drugs = drugs, universe = "Default", 
##' #                                 type = "GO", ont="MF")
##' # dup_hyperG_res 

setMethod("show", signature(object="feaResult"),
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

##' @name result
##' @docType methods
##' @rdname result-methods
##' @method result feaResult
##' @aliases result,feaResult-method
##' @examples
##' data(drugs)
##' #dup_hyperG_res <- tsea_dup_hyperG(drugs = drugs, universe = "Default", 
##' #                                 type = "GO", ont="MF")
##' #result(dup_hyperG_res) 

setMethod("result", signature(x="feaResult"),
          function(x) x@result)

##' @description The \code{get_drugs} generic extracts the drug names/ids stored
##' in the \code{drugs} slot of an feaResult object.
##' @name get_drugs
##' @docType methods
##' @rdname get_drugs-methods
##' @method get_drugs feaResult
##' @aliases get_drugs,feaResult-method
##' @aliases get_drugs-method
##' @param x feaResult object
##' @return character vector
##' @examples 
##' data(drugs)
##' # dup_hyperG_res <- tsea_dup_hyperG(drugs = drugs, universe = "Default", 
##' #                                   type = "GO", ont="MF")
##' # get_drugs(dup_hyperG_res)

setMethod("get_drugs", signature(x="feaResult"),
          function(x) x@drugs)

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
# 
# ##' cnetplot for feaResult
# ##' 
# ##' @rdname cnetplot-methods
# ##' @aliases cnetplot,feaResult,ANY-method
# ##' @param x feaResult object
# ##' @param showCategory number of enriched terms to display
# ##' @param categorySize one of 'pvalue', 'p.adjust' and 'qvalue'
# ##' @param foldChange fold change
# ##' @param fixed TRUE or FALSE
# ##' @importFrom enrichplot cnetplot
# ##' @export
# setMethod("cnetplot", signature(x="feaResult"),
#           function(x, showCategory=5, categorySize="pvalue", 
#                    foldChange=NULL, fixed=TRUE, ...) {
#               cnetplot.feaResult(x,
#                                     showCategory=showCategory,
#                                     categorySize=categorySize,
#                                     foldChange=foldChange,
#                                     fixed=fixed, ...)
# })

