#' show \code{\link{qSig}}, \code{\link{gessResult}}, \code{\link{feaResult}}
#' objects
#' 
#' @name show
#' @docType methods
#' @rdname show-methods
#' @title show method
#' @param object object used for show
#' @return message
#' @aliases show,gessResult-method
#' @usage show(object)
#' @examples 
#' gr <- gessResult(result=dplyr::tibble(pert=letters[seq_len(10)], 
#'                                val=seq_len(10)), 
#'                  query=list(up=c("g1","g2"), down=c("g3","g4")),
#'                  gess_method="LINCS", refdb="path/to/lincs/db")
#' gr
setMethod("show", signature(object="gessResult"),
          function (object) {
            cat("#\n# gessResult object \n#\n")
            cat("@result \n")
            print(result(object))
            q <- qr(object)
            if(is(q, "list")){
                if(length(q$upset)>10){
                   cat("@query", "\t", "up gene set", 
                        paste0("(", length(q$upset), "):"), 
                        "\t", q$upset[seq_len(10)], "... \n")
                 } else {
                     cat("@query", "\t", "up gene set", 
                         paste0("(", length(q$upset), "):"), 
                         "\t", q$upset, "\n")
                }
                if(length(q$downset)>10){
                    cat("     ", "\t", "down gene set", 
                        paste0("(", length(q$downset), "):"), 
                        "\t", q$downset[seq_len(10)], "... \n")
                } else {
                    cat("     ", "\t", "down gene set", 
                        paste0("(", length(q$downset), "):"), 
                        "\t", q$downset, "\n")
                }
            }
            if(is(q, "matrix")){
                cat("@query\n")
                if(nrow(q)>10){
                    print(head(q,10))
                    cat("# ... with", nrow(q)-10, "more rows\n")
                } else {
                    print(q)
                }
            }
            cat("\n@gess_method", "\t", gm(object), "\n")
            cat("\n@refdb", "\t", refdb(object), "\n")
          })

#' @name result
#' @docType methods
#' @rdname result-methods
#' @method result gessResult
#' @param x \code{gessResult} or \code{feaResult} object
#' @return tibble
#' @aliases result,gessResult-method
#' @examples 
#' gr <- gessResult(result=dplyr::tibble(pert=letters[seq_len(10)], 
#'                                val=seq_len(10)), 
#'                  query=list(up=c("g1","g2"), down=c("g3","g4")),
#'                  gess_method="LINCS", refdb="path/to/lincs/db")
#' result(gr)

setMethod("result", c(x="gessResult"),
          function(x) x@result)

