##'
##' It builds a `qSig` object to store the query signature, reference database
##' and GESS method used for GESS methods
##' @docType methods
##' @name qSig
##' @rdname qSig-methods
##' @aliases qSig,list,character,character-method
##' @param query If 'gess_method' is 'CMAP' or 'LINCS', it should be a list with
##' two character vectors named \code{upset} and \code{downset} for up- and
##' down-regulated gene labels, respectively. The labels should be gene Entrez
##' IDs if the reference database is a pre-built CMAP or LINCS database. If a
##' custom database is used, the labels need to be of the same type as those in
##' the reference database.
##' 
##' If 'gess_method' is 'gCMAP', the query is a matrix with a single column
##' representing gene ranks from a biological state of interest. The corresponding 
##' gene labels are stored in the row name slot of the matrix. Instead of ranks
##' one can provide scores (e.g. z-scores). In such a case, the scores will be 
##' internally transformed to ranks. 
##' 
##' If 'gess_method' is 'Fisher', the query is expected to be a list with two
##' character vectors named \code{upset} and \code{downset} for up- and
##' down-regulated gene labels, respectively (same as for 'CMAP' or 'LINCS' method).
##' Internally, the up/down gene labels are combined into a single gene set
##' when querying the reference database with the Fisher's exact test. This means the
##' query is performed with an unsigned set. The query can also be a matrix with a
##' single numeric column and the gene labels (e.g. Entrez gene IDs) in the row name
##' slot. The values in this matrix can be z-scores or LFCs. In this case, the
##' actual query gene set is obtained according to upper and lower cutoffs set by
##' the user.
##' 
##' If 'gess_method' is 'Cor', the query is a matrix with a single numeric column
##' and the gene labels in the row name slot. The numeric column can contain z-scores,
##' LFCs, (normalized) gene expression intensity values or read counts.
##' @param gess_method one of 'CMAP', 'LINCS', 'gCMAP', 'Fisher' or 'Cor'
##' @param refdb character(1), can be one of "cmap", "cmap_expr", "lincs", or 
##' "lincs_expr" when using the CMAP/LINCS databases from the affiliated
##' \code{signatureSearchData} package. With 'cmap' the database contains
##' signatures of LFC scores obtained from DEG analysis routines; with
##' 'cmap_expr' normalized gene expression values; with 'lincs' z-scores
##' obtained from the DEG analysis methods of the LINCS project; and with
##' 'lincs_expr' normalized expression values.
##' 
##' To use a custom signature database, it should be the file path to the HDF5 file
##' generated with the \code{\link{build_custom_db}} function.  Alternatively, a
##' suitable version of the CMAP/LINCS databases can be used. For details on this,
##' please consult the vignette of the \pkg{signatureSearchData} package.
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
##' qsig_lincs <- qSig(query=list(upset=upset, downset=downset), 
##'                    gess_method="LINCS", refdb=db_path)
##' qsig_gcmap <- qSig(query=query_mat, gess_method="gCMAP", refdb=db_path)
##' @exportMethod qSig
setMethod("qSig",
  signature(query="list", gess_method="character", 
            refdb="character"),
  function(query, gess_method, refdb){
    ## Validity check of refdb
    if(!any(refdb %in% c("cmap","cmap_expr","lincs","lincs_expr"))){
        ref_val <- h5read(refdb, "assay", c(1,1))
        if(!is.numeric(ref_val)) 
            stop("The value stored in 'refdb' should be numeric!")
    }
    if(any(gess_method %in% c("CMAP", "LINCS", "Fisher"))){
      upset = query[[1]]
      downset = query[[2]]
      refdb = determine_refdb(refdb)
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
##' @aliases qSig,matrix,character,character-method
##' @exportMethod qSig
setMethod("qSig",
  signature(query="matrix", gess_method="character", refdb="character"),
  function(query, gess_method, refdb){
    ## Validity check of refdb
    if(!any(refdb %in% c("cmap","cmap_expr","lincs","lincs_expr"))){
      ref_val <- h5read(refdb, "assay", c(1,1))
      if(!is.numeric(ref_val)) 
        stop("The value stored in 'refdb' should be numeric!")
    }
    if(any(gess_method %in% c("gCMAP", "Fisher", "Cor"))){
      experiment = query
      refdb = determine_refdb(refdb)
      gid_db <- h5read(refdb, "rownames", drop=TRUE)
      if(! is.numeric(experiment[1,1])) 
        stop("The 'qsig' should be a numeric matrix 
             representing genome-wide signature from a treatment")
      if(is.null(rownames(experiment))) 
        stop("The gene labels should be stored in the row name slot of the 
             matrix!")
      if(sum(rownames(experiment) %in% gid_db)==0) 
        stop("The gene labels in the query matrix share 0 identifiers with 
             reference database!")
    } else
      stop("'gess_method' slot must be one of 'gCMAP', 'Fisher' or 'Cor' 
  if 'qsig' is a numeric matrix with gene labels in the row name slot!")
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

