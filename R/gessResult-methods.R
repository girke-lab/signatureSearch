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
##' @examples 
##' db_path <- system.file("extdata", "sample_db.h5", 
##' package = "signatureSearch")
##' # Load sample_db as `SummarizedExperiment` object
##' library(signatureSearchData)
##' sample_db <- readHDF5chunk(db_path, colindex=1:100)
##' # get "vorinostat__SKB__trt_cp" signature drawn from sample databass
##' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
##' qsig_fisher <- qSig(query=query_mat, gess_method="Fisher", refdb=db_path)
##' fisher <- gess_fisher(qSig=qsig_fisher, higher=1, lower=-1)
##' fisher
setMethod("show", signature(object="gessResult"),
          function (object) {
            cat("#\n# gessResult object \n#\n")
            cat("@result \n")
            print(object@result)
            if(is(object@query, "list")){
              cat("@query", "\t", "up gene set", 
                  paste0("(", length(object@query[[1]]), "):"), 
                  "\t", object@query[[1]][seq_len(10)], "... \n")
              cat("     ", "\t", "down gene set", 
                  paste0("(", length(object@query[[2]]), "):"), 
                  "\t", object@query[[2]][seq_len(10)], "... \n")
            }
            if(is(object@query, "matrix")){
              cat("@query\n")
              mat=object@query
              print(head(mat,10))
              cat("# ... with", nrow(mat)-10, "more rows\n")
            }
            cat("\n@gess_method", "\t", object@gess_method, "\n")
            cat("\n@refdb", "\t", object@refdb, "\n")
          })

##' @description get 'result' slot of gessResult object
##' @name result
##' @docType methods
##' @rdname result-methods
##' @method result gessResult
##' @param x \code{gessResult} or \code{feaResult} object
##' @return tibble
##' @aliases result,gessResult-method
##' @examples 
##' db_path <- system.file("extdata", "sample_db.h5", 
##' package = "signatureSearch")
##' # Load sample_db as `SummarizedExperiment` object
##' library(signatureSearchData)
##' sample_db <- readHDF5chunk(db_path, colindex=1:100)
##' # get "vorinostat__SKB__trt_cp" signature drawn from sample databass
##' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
##' qsig_fisher <- qSig(query=query_mat, gess_method="Fisher", refdb=db_path)
##' fisher <- gess_fisher(qSig=qsig_fisher, higher=1, lower=-1)
##' result(fisher)

setMethod("result", signature(x="gessResult"),
          function(x) x@result)

