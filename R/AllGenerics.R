##' cnetplot
##'
##' @docType methods
##' @name cnetplot
##' @rdname cnetplot-methods
##' @title cnetplot method
##' @param x enrichResult object
##' @param showCategory number of category plotted
##' @param categorySize one of geneNum or pvalue
##' @param foldChange fold change of expression value
##' @param fixed logical
##' @param ... additional parameters
##' @return plot
##' @export
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
setGeneric("cnetplot",
           function(x, showCategory=5, categorySize="geneNum", foldChange=NULL, fixed=TRUE, ...)
               standardGeneric("cnetplot"))

##' dotplot
##'
##' @docType methods
##' @name dotplot
##' @rdname dotplot-methods
##' @title dotplot method
##' @param ... additional parameter
##' @return plot
##' @export
##' @author Guangchuang Yu
setGeneric("dotplot", function(object, ...) standardGeneric("dotplot"))

#' geneID generic
#'
#' @param x feaResult object
#' @return 'geneID' return the 'geneID' column of the FEA result which can be converted to data.frame via 'as.data.frame'
#' @export
#' @examples
#' \dontrun{
#' data(geneList, package="DOSE")
#' de <- names(geneList)[1:100]
#' x <- enrichDO(de)
#' geneID(x)
#' }

geneID <- function(x) {
   UseMethod("geneID", x)
}

#' geneInCategory generic
#'
#' @param x feaResult
#' @return 'geneInCategory' return a list of genes, by spliting the input gene vector to enriched functional categories
#' @export
#' @examples
#' \dontrun{
#' data(geneList, package="DOSE")
#' de <- names(geneList)[1:100]
#' x <- enrichDO(de)
#' geneInCategory(x)
#' }

geneInCategory <- function(x) {
   UseMethod("geneInCategory", x)
}

##' qSig generic
##' 
##' @title qSig method
##' @rdname qSig-methods
##' @export
setGeneric("qSig", function(qsig, gess_method, refdb) standardGeneric("qSig"))

##' result generic
##' 
##' @title result method
##' @rdname result-methods
##' @export
setGeneric("result", function(x) standardGeneric("result"))



