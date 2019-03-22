##' qSig method
##'
##' It builds a `qSig` object to store the query signature, reference database
##' and GESS method used to search for similarity
##' @docType methods
##' @name qSig
##' @rdname qSig-methods
##' @aliases qSig,list,character,SummarizedExperiment-method
##' @param qsig When 'gess_method' is 'CMAP' or 'LINCS', it should be a list of 
##' two elements, which are up and down regulated gene sets of entrez ids.
##' 
##' When 'gess_method' is 'gCMAP', 'Fisher' or 'Cor', it should be a matrix 
##' representing gene expression profiles (GEPs) of treatment(s). 
##' @param gess_method one of 'CMAP', 'LINCS', 'gCMAP', 'Fisher' or 'Cor'
##' @param refdb \code{SummarizedExperiment} object, which can be HDF5 backed 
##' and loaded via `loadHDF5SummarizedExperiment` function. 
##' The 'assays' slot of the \strong{SummarizedExperiment} object should be a 
##' \code{DelayedMatrix} or a matrix consists of genome-wide (GEPs) from a
##' number of drug treatments or genetic perturbations. It represents the 
##' reference database that the query signature is searched against. 
##' 
##' The \code{sample_db} contains 95 GEPs randomly sampled from the `lincs` 
##' database and 5 GEPs from HDAC inhibitors in human SKB (muscle) cell. 
##' 
##' The full `lincs` and `cmap` public databases can be loaded from the 
##' \code{\link{signatureSerch_data}} package.
##' 
##' The custom database can be built via \code{\link{`build_custom_db`}} 
##' function if a `data.frame` representing genome-wide GEPs (log2FC, z-scores, 
##' intensity values, etc.) of compound or genetic treatments in cells 
##' is provided.
##' @return \code{qSig} object
##' @seealso \code{\link{build_custom_db}}, \code{\link{signatureSerch_data}},
##'          \code{\link[SummarizedExperiment]{SummarizedExperiment}},
##'          \code{\link[DelayedArray]{DelayedMatrix}},
##'          \code{\link[HDF5Array]{loadHDF5SummarizedExperiment}}
##' @examples 
##' db_dir <- system.file("extdata", "sample_db", package = "signatureSearch")
##' sample_db <- loadHDF5SummarizedExperiment(db_dir)
##' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
##' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
##' query = as.numeric(query_mat); names(query) = rownames(query_mat)
##' upset <- head(names(query[order(-query)]), 150)
##' downset <- tail(names(query[order(-query)]), 150)
##' qsig_lincs <- qSig(qsig = list(upset=upset, downset=downset), 
##'                    gess_method = "LINCS", refdb = sample_db)
##' qsig_gcmap <- qSig(qsig=query_mat, gess_method="gCMAP", refdb=sample_db)
##' @exportMethod qSig
setMethod("qSig",
  signature(qsig="list", gess_method="character", refdb="SummarizedExperiment"),
  function(qsig, gess_method, refdb){
    ## Validity check of refdb
    if(!is.numeric(as.matrix(assay(refdb)[1,1]))) 
      stop("The value stored in 'assays' slot of 'refdb' should be numeric!")
    if(any(gess_method %in% c("CMAP", "LINCS"))){
      upset = qsig[[1]]
      downset = qsig[[2]]
      se = refdb
      gid_db <- rownames(se)
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
      qsig[[1]] = upset
      qsig[[2]] = downset
    } else
      stop("'gess_method' slot must be one of 'CMAP', 'LINCS', or 'Fisher' if 
'qsig' is a list of two elements representing up and down regulated gene sets!")
    new("qSig", qsig=qsig, gess_method=gess_method, refdb=refdb)
  }
)

##' qSig method
##'
##' @docType methods
##' @name qSig
##' @rdname qSig-methods
##' @aliases qSig,matrix,character,SummarizedExperiment-method
##' @exportMethod qSig
setMethod("qSig",
  signature(qsig="matrix", gess_method="character", 
            refdb="SummarizedExperiment"),
  function(qsig, gess_method, refdb){
    ## Validity check of refdb
    if(!is.numeric(as.matrix(assay(refdb)[1,1]))) 
      stop("The value stored in 'assays' slot of 'refdb' should be numeric")
    if(any(gess_method %in% c("gCMAP", "Fisher", "Cor"))){
      experiment = qsig
      gid_db <- rownames(refdb)
      if(! is.numeric(experiment[1,1])) 
        stop("The 'qsig' should be a numeric matrix 
             representing genome-wide GEPs from treatments")
      if(is.null(rownames(experiment))) 
        stop("The 'qsig' should be a numeric matrix with rownames as gene 
             Entrez ids representing genome-wide GEPs")
      if(sum(rownames(experiment) %in% gid_db)==0) 
        stop("The rownames of 'qsig' share 0 identifiers with 
             reference database!")
    } else
      stop("'gess_method' slot must be one of 'gCMAP', 'Fisher' or 'Cor' 
  if 'qsig' is a numeric matrix representing genome-wide GEPs from treatments!")
    new("qSig", qsig=qsig, gess_method=gess_method, refdb=refdb)
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

