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
              # cat("#...Citation\n")
              # citation_msg <- paste("  Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.",
              #                       "  clusterProfiler: an R package for comparing biological themes among",
              #                       "  gene clusters. OMICS: A Journal of Integrative Biology",
              #                       "  2012, 16(5):284-287", sep="\n", collapse="\n")
              # cat(citation_msg, "\n\n")
          })

##' plot method generics
##'
##' @docType methods
##' @name plot
##' @rdname plot-methods
##' @aliases plot,feaResult,ANY-method
##' @title plot method
##' @param x A \code{feaResult} instance
##' @param type one of dot, bar, cnet or enrichMap
##' @param ... Additional argument list
##' @return plot
##' @importFrom stats4 plot
##' @exportMethod plot
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
setMethod("plot", signature(x="feaResult"),
          function(x, type = "bar", ... ) {
              if (type == "cnet" || type == "cnetplot") {
                  cnetplot.enrichResult(x, ...)
              }
              if (type == "bar" || type == "barplot") {
                  barplot(x, ...)
              }
              if (type == "enrichMap") {
                  enrichMap(x, ...)
              }
              if (type == "dot" || type == "dotplot") {
                  dotplot(x, ...)
              }
          }
          )


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
##' @exportMethod dotplot
##' @author Guangchuang Yu
setMethod("dotplot", signature(object="feaResult"),
          function(object, x="geneRatio", colorBy="p.adjust", showCategory=10, split=NULL, font.size=12, title="") {
              dotplot_internal(object, x, colorBy, showCategory, split, font.size, title)
          }
          )

##' @rdname cnetplot-methods
##' @exportMethod cnetplot
setMethod("cnetplot", signature(x="feaResult"),
          function(x, showCategory=5, categorySize="pvalue", foldChange=NULL, fixed=TRUE, ...) {
              cnetplot.enrichResult(x,
                                    showCategory=showCategory,
                                    categorySize=categorySize,
                                    foldChange=foldChange,
                                    fixed=fixed, ...)
          }
          )
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

