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
              cat("@qsig", "\t", "up gene set", 
                  paste0("(", length(object@qsig[[1]]), "):"), 
                  "\t", object@qsig[[1]][seq_len(10)], "... \n")
              cat("     ", "\t", "down gene set", 
                  paste0("(", length(object@qsig[[2]]), "):"), 
                  "\t", object@qsig[[2]][seq_len(10)], "... \n")
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

##' @name show
##' @docType methods
##' @rdname show-methods
##' @aliases show,feaResult-method
##' @importFrom utils str
setMethod("show", signature(object="feaResult"),
          function (object){
              cat("#\n# Functional Enrichment Analysis \n#\n")
              cat("#...@organism", "\t", object@organism, "\n")
              cat("#...@ontology", "\t", object@ontology, "\n")
              cat("#...@drugs", "\t")
              str(object@drugs)
              cat("#...@targets", "\t")
              str(object@targets)
              cat("#...@universe", "\t")
              str(object@universe)
              cat(paste0("#...", nrow(object@result)), "enriched terms found\n")
              str(object@result)
          })


##' @description get 'result' slot of gessResult object
##' @name result
##' @docType methods
##' @rdname result-methods
##' @method result gessResult
##' @param x \code{gessResult} or \code{feaResult} object
##' @return tibble
##' @aliases result,gessResult-method
setMethod("result", signature(x="gessResult"),
          function(x) x@result)


##' @description get 'result' slot of feaResult object
##' @name result
##' @docType methods
##' @rdname result-methods
##' @method result feaResult
##' @aliases result,feaResult-method
setMethod("result", signature(x="feaResult"),
          function(x) x@result)

##' @description get `drugs` slot of feaResult object
##' @name get_drugs
##' @docType methods
##' @rdname get_drugs-methods
##' @method get_drugs feaResult
##' @aliases get_drugs,feaResult-method
##' @aliases get_drugs-method
##' @param x feaResult object
##' @return character vector
setMethod("get_drugs", signature(x="feaResult"),
          function(x) x@drugs)

##' dotplot for feaResult
##'
##' @rdname dotplot-methods
##' @aliases dotplot,feaResult,ANY-method
##' @param object an instance of feaResult
##' @param x variable for x axis
##' @param colorBy one of 'pvalue', 'p.adjust' and 'qvalue'
##' @param showCategory number of category
##' @param split separate result by 'category' variable
##' @param font.size font size
##' @param title plot title
##' @importFrom enrichplot dotplot
##' @exportMethod dotplot
setMethod("dotplot", signature(object="feaResult"),
          function(object, x="geneRatio", colorBy="p.adjust", showCategory=10, 
                   split=NULL, font.size=12, title="") {
              dotplot_internal(object, x, colorBy, showCategory, 
                               split, font.size, title)
          })

##' cnetplot for feaResult
##' 
##' @rdname cnetplot-methods
##' @aliases cnetplot,feaResult,ANY-method
##' @param x feaResult object
##' @param showCategory number of enriched terms to display
##' @param categorySize one of 'pvalue', 'p.adjust' and 'qvalue'
##' @param foldChange fold change
##' @param fixed TRUE or FALSE
##' @importFrom enrichplot cnetplot
##' @export
setMethod("cnetplot", signature(x="feaResult"),
          function(x, showCategory=5, categorySize="pvalue", 
                   foldChange=NULL, fixed=TRUE, ...) {
              cnetplot.feaResult(x,
                                    showCategory=showCategory,
                                    categorySize=categorySize,
                                    foldChange=foldChange,
                                    fixed=fixed, ...)
})
# 
# ##' @rdname dtnetplot-methods
# ##' @exportMethod dtnetplot
# setMethod("dtnetplot", signature(x="enrichResult"),
#           function(c_ego, GOterm, drugSize="targetNum", targetWeight=NULL, fixed=TRUE, ...) {
#             dtnetplot.enrichResult(c_ego,
#                                    GOterm,
#                                    drugSize=drugSize,
#                                    targetWeight=targetWeight,
#                                    fixed=fixed, ...)
#           }
# )

