##' qSig method
##'
##' It builds a `qSig` object to store the query signature, reference database
##' and GESS method used to search for similarity
##' @docType methods
##' @name qSig
##' @rdname qSig-methods
##' @aliases qSig,list,character,character,character-method
##' @param qsig When 'gess_method' is 'CMAP' or 'LINCS', 
##' it should be a list of two elements, which are up and down regulated gene 
##' sets of entrez ids.
##' 
##' When 'gess_method' is 'gCMAP', 'Fisher' or 'Cor', it should be a matrix 
##' representing gene expression profiles (GEPs) of treatment(s). 
##' @param gess_method one of 'CMAP', 'LINCS', 'gCMAP', 'Fisher' or 'Cor'
##' @param refdb character(1), can be "cmap", "cmap_expr", "lincs", or 
##' "lincs_expr" if users want to use the existing CMAP/LINCS databases. 
##' 
##' If users want to use the custom signature database, 
##' it should be the file path to the HDF5 file generated with 
##' \code{\link{build_custom_db}} function or
##' generated from the source files of CMAP/LINCS databases according to 
##' the vignette in \pkg{signatureSearchData} package. The HDF5 file contains 
##' the reference signatures that the query signature is searched against. 
##' @return \code{qSig} object
##' @importFrom rhdf5 h5read
##' @seealso \code{\link{build_custom_db}}, 
##' \code{\link[signatureSearchData]{signatureSearchData}}
##' @examples 
##' db_path <- system.file("extdata", "sample_db.h5", 
##' package = "signatureSearch")
##' # Load sample_db as `SummarizedExperiment` object
##' library(signatureSearchData)
##' sample_db <- readHDF5chunk(db_path, colindex=1:100)
##' # get "vorinostat__SKB__trt_cp" signature drawn from sample database
##' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
##' query = as.numeric(query_mat); names(query) = rownames(query_mat)
##' upset <- head(names(query[order(-query)]), 150)
##' downset <- tail(names(query[order(-query)]), 150)
##' qsig_lincs <- qSig(query = list(upset=upset, downset=downset), 
##'                    gess_method="LINCS", refdb=db_path)
##' qsig_gcmap <- qSig(qsig=query_mat, gess_method="gCMAP", refdb=db_path)
##' @exportMethod qSig
setMethod("qSig",
  signature(query="list", gess_method="character", 
            refdb="character"),
  function(query, gess_method, refdb){
    ## Validity check of refdb
    ref_val <- h5read(refdb, "assay", c(1,1))
    if(!is.numeric(ref_val)) 
      stop("The value stored in 'refdb' should be numeric!")
    if(any(gess_method %in% c("CMAP", "LINCS"))){
      upset = query[[1]]
      downset = query[[2]]
      gid_db <- h5read(refdb, "rownames", drop=TRUE)
      ## Validity checks of upset and downset
      if(all(c(!is.character(upset), !is.null(upset)))) 
        stop("upset of 'qsig' slot needs to be ID character vector or NULL")
      if(all(c(!is.character(downset), !is.null(downset)))) 
        stop("downset of 'qsig' slot needs to be ID character vector or NULL")
      if(is.null(upset) & is.null(downset)) 
        stop("both or one of the upset and downset in 'qsig' slot need to be 
             assigned query entrez IDs as character vector")
      
      ## Remove entries in up/down set not present in reference database
      if(!is.null(upset)){
        message(paste(sum(upset %in% gid_db), "/", length(upset), 
                  "genes in up set share identifiers with reference database"))
        upset <- upset[upset %in% gid_db]
        if(length(upset)==0) 
          stop("upset shares zero idenifiers with reference database, 
               please set upset of 'qsig' slot as NULL")
      }
      if(!is.null(downset)){
        message(paste(sum(downset %in% gid_db),"/",length(downset), 
                "genes in down set share identifiers with reference database"))
        downset <- downset[downset %in% gid_db]
        if(length(downset)==0) 
          stop("downset shares zero idenifiers with reference database, 
               please set downset of 'qsig' slot as NULL")
      }
      query[[1]] = upset
      query[[2]] = downset
    } else
      stop("'gess_method' slot must be one of 'CMAP', 'LINCS', or 'Fisher' if 
'qsig' is a list of two elements representing up and down regulated gene sets!")
    new("qSig", query=query, gess_method=gess_method, refdb=refdb)
  }
)

##' qSig method
##'
##' @docType methods
##' @name qSig
##' @rdname qSig-methods
##' @aliases qSig,matrix,character,character,character-method
##' @exportMethod qSig
setMethod("qSig",
  signature(query="matrix", gess_method="character", refdb="character"),
  function(query, gess_method, refdb){
    ## Validity check of refdb
    ref_val <- h5read(db_path, "assay", c(1,1))
    if(!is.numeric(ref_val)) 
      stop("The value stored in 'assays' slot of 'refdb' should be numeric")
    if(any(gess_method %in% c("gCMAP", "Fisher", "Cor"))){
      experiment = query
      gid_db <- h5read(refdb, "rownames", drop=TRUE)
      if(! is.numeric(experiment[1,1])) 
        stop("The 'qsig' should be a numeric matrix 
             representing genome-wide GEPs from treatments")
      if(is.null(rownames(experiment))) 
        stop("The 'query' should be a numeric matrix with rownames as gene 
             Entrez ids representing genome-wide GEPs")
      if(sum(rownames(experiment) %in% gid_db)==0) 
        stop("The rownames of 'query' share 0 identifiers with 
             reference database!")
    } else
      stop("'gess_method' slot must be one of 'gCMAP', 'Fisher' or 'Cor' 
  if 'qsig' is a numeric matrix representing genome-wide GEPs of treatments!")
    new("qSig", query=query, gess_method=gess_method, refdb=refdb)
  }
)

##' @name show
##' @docType methods
##' @rdname show-methods
##' @aliases show,qSig-method
##' @exportMethod show
setMethod("show", signature(object="qSig"),
    function (object) {
      cat("#\n# qSig object used for GESS analysis \n#\n")
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
      cat("\n@refdb", "\t")
      cat(object@refdb, "\n")
    })

