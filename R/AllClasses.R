##' Class "qSig"
##' 
##' This class stores the query signature, reference database
##' and GESS method used to search for similarity
##' @name qSig-class
##' @docType class
##' @aliases qSig-class
##' @slot qsig When 'gess_method' is 'CMAP' or 'LINCS', it should be a list of 
##' two elements, which are up and down regulated gene sets of entrez ids.
##' 
##' When 'gess_method' is 'gCMAP', 'Fisher' or 'Cor', it should be a matrix 
##' representing gene expression profiles (GEPs) of treatment(s). 
##' @slot gess_method one of 'CMAP', 'LINCS', 'gCMAP', 'Fisher' or 'Cor'
##' @slot refdb \code{SummarizedExperiment} object, which can be HDF5 backed 
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
##' @slot refdb_name character, name of the reference database. Like "CMAP",
##' "LINCS" or other custom names.
##' @exportClass qSig
##' @keywords classes
setClass("qSig", slots = c(
  qsig = "ANY",
  gess_method = "character",
  refdb = "SummarizedExperiment",
  refdb_name = "character"
))

##' Class "gessResult"
##' 
##' The class stores the result of GESS analysis
##' @name gessResult-class
##' @aliases gessResult
##' @docType class
##' @slot result tibble from GESS analysis, represents a list of drugs 
##' in the reference database ranked by their signature similarity to the query
##' @slot qsig query signature
##' @slot gess_method method for GESS analysis
##' @slot refdb_name name of the reference database
##' @exportClass gessResult
##' @keywords classes
setClass("gessResult",
         slots = c(
           result = "data.frame",
           qsig = "ANY",
           gess_method = "character",
           refdb_name = "character"
         ))

## Constructor for "gessResult"
gessResult <- function(result, qsig, gess_method, refdb_name="UNKNOWN")
  new("gessResult", result=result, qsig=qsig, 
      gess_method=gess_method, refdb_name=refdb_name)


## Defining the validity method for "qSig"
# setValidity("qSig", function(object) {
#      TRUE
# })

##' Class "feaResult"
##' 
##' This class represents the result of functional enrichment analysis.
##'
##'
##' @name feaResult-class
##' @docType class
##' @slot result data.frame from FEA analysis
##' @slot organism only "human" supported
##' @slot ontology biological ontology
##' @slot drugs Drug IDs
##' @slot targets Target IDs of drugs in DrugBank/LINCS/STITCH databases or 
##' target list with scores.
##' @exportClass feaResult
##' @author Yuzhu Duan
##' @keywords classes
setClass("feaResult",
         representation=representation(
           result         = "data.frame",
           organism       = "character",
           ontology       = "character",
           drugs          = "character",
           targets        = "ANY"
         )
)
## @slot universe background genes or drugs. For TSEA, it is all the genes 
## in the corresponding annotation system (GO/KEGG). For DSEA, it is all the 
## drugs in the correspoinding annotation system (GO/KEGG) after 
## drug-to-functional category mapping