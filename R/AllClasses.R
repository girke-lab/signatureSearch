##' Class "qSig"
##' 
##' This class represents the query signature for GESS analysis
##' @name qSig-class
##' @docType class
##' @slot qsig When "gess_method" is "CMAP", "LINCS" or "Fisher", "sig" is a list of two elements, up and down regulated gene sets of entrez ids. 
##' When "gess_method" is "gCMAP" or "SP", "sig" is a list of one element, named numeric vector represents GEP from DE analysis or gene expression values.
##' @slot gess_method one of "CMAP","LINCS","gCMAP","Fisher","SP"
##' @slot refdb \code{SummarizedExperiment} object or a HDF5 backend \code{SummarizedExperiment} object loaded via 
##' `loadHDF5SummarizedExperiment` function from the HDF5 datasets saved in a directory. 
##' where all the assays are \code{HDF5Array} or \code{DelayedArray} objects consist of genome-wide differential expression profiles (GEPs)
##'  (e.g. log2 ratios or z-scores) from various drug treatments or genetic perturbations. 
##' It represents reference database that the query signature is searched against. Can be existing public databases (CMAP or LINCS)
##' or custom database.
##' @exportClass qSig
##' @keywords classes
setClass("qSig", slots = c(
  qsig = "ANY",
  gess_method = "character",
  refdb = "SummarizedExperiment"
))
  
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
##' @slot drug Drug IDs
##' @slot universe background genes or drugs. For TSEA, it is all the genes in the corresponding annotation system (GO/KEGG). For DSEA, it is all the drugs
##' in the correspoinding annotation system (GO/KEGG) after drug-to-functional category mapping
##' @slot refSets gene sets or drug sets in the corresponding annotation system
##' @exportClass feaResult
##' @author Yuzhu Duan
##' @keywords classes
setClass("feaResult",
         representation=representation(
             result         = "data.frame",
             organism       = "character",
             ontology       = "character",
             drug           = "character",
             universe       = "character",
             refSets        = "list"
             )
         )

##' Class "gessResult"
##' 
##' The class represents the result of GESS analysis
##' @name gessResult-class
##' @docType class
##' @slot result tibble from GESS analysis
##' @slot qsig query signature
##' @slot gess_method method for GESS analysis
##' @slot refdb `SummarizedExperiment` object represents refrence signaure database
##' @exportClass gessResult
##' @keywords classes
setClass("gessResult",
         slots = c(
           result = "data.frame",
           qsig = "ANY",
           gess_method = "character",
           refdb = "SummarizedExperiment"
         ))
