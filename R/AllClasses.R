##' Class "qSig"
##' 
##' This class stores the query signature, reference database
##' and GESS method used to search for similarity
##' @name qSig-class
##' @docType class
##' @aliases qSig-class
##' @slot query When 'gess_method' is 'CMAP' or 'LINCS', 
##' it should be a list of two elements, which are up and down regulated 
##' gene sets of entrez ids.
##' 
##' When 'gess_method' is 'gCMAP', 'Fisher' or 'Cor', it should be a matrix 
##' representing gene expression profiles (GEPs) of treatment(s). 
##' @slot gess_method one of 'CMAP', 'LINCS', 'gCMAP', 'Fisher' or 'Cor'
##' @slot refdb character(1), can be "cmap", "cmap_expr", "lincs", or 
##' "lincs_expr" if users want to use the existing CMAP/LINCS databases. 
##' 
##' If users want to use the custom signature database, 
##' it should be the file path to the HDF5 file generated with 
##' \code{\link{build_custom_db}} function or
##' generated from the source files of CMAP/LINCS databases according to 
##' the vignette in \code{\link[signatureSearchData]{signatureSearchData}}
##' package. The HDF5 file contains 
##' the reference signatures that the query signature is searched against. 
##' @exportClass qSig
##' @keywords classes
setClass("qSig", slots = c(
  query = "ANY",
  gess_method = "character",
  refdb = "character"
))

##' gessResult object
##' 
##' The gessResult object stores the search result table, query signature, 
##' name of the GESS method and path to the reference database from the 
##' GESS methods.
##' @name gessResult-class
##' @aliases gessResult
##' @docType class
##' @slot result tibble object, this result table contains the search results 
##' for each perturbagen in the reference database ranked by their signature 
##' similarity to the query. The result table can be extracted via 
##' \code{\link{result}} accessor function.
##' 
##' Description of the common columns from different GESS methods:
##' \itemize{
##'     \item pert: character, name of perturbagen (e.g. drug) in the reference 
##'     database
##'     \item cell: character, acronym of cell type
##'     \item type: character, perturbation type. In CMAP and LINCS 
##'     databases, the perturbation types are all compound treatment(trt_cp). 
##'     Users can build their custom signature database with other types of 
##'     perturbation, e.g., gene knockdown or overexpression via 
##'     \code{\link{build_custom_db}} function
##'     \item trend: character, up or down when reference signature is 
##'     positively or negatively connected with the query signature, 
##'     respectively.
##'     \item N_upset: integer, number of genes in the query up set
##'     \item N_downset: integer, number of genes in the query down set
##'     \item t_gn_sym: character, gene symbols of the corresponding drug 
##'     targets
##' } 
##' @slot query query signature
##' @slot gess_method name of the GESS method 
##' @slot refdb path to the reference database
##' @exportClass gessResult
##' @keywords classes
setClass("gessResult",
         slots = c(
           result = "data.frame",
           query = "ANY",
           gess_method = "character",
           refdb = "character"
         ))

## Constructor for "gessResult"
gessResult <- function(result, query, gess_method, refdb)
  new("gessResult", result=result, query=query, 
      gess_method=gess_method, refdb=refdb)


## Defining the validity method for "qSig"
# setValidity("qSig", function(object) {
#      TRUE
# })

##' feaResult object
##' 
##' The feaResult object stores the enrichment result table, organism 
##' information of the annotation system, and the ontology type of the GO 
##' annotation system. If the annotation system is KEGG, the latter will be 
##' 'KEGG'. It also stores the input drugs used for the enrichment test, 
##' as well as their target information.
##' @name feaResult-class
##' @aliases feaResult
##' @docType class
##' @slot result tibble object, this result table contains the
##' enriched functional categories (e.g. GO terms or KEGG pathways) ranked by 
##' the corresponding enrichment statistic. The result table can be extracted 
##' via \code{\link{result}} accessor function.
##' 
##' Description of the common columns in the result table from different 
##' enrichment methods:
##' \itemize{
##'     \item ont: in case of GO, one of BP, MF, CC, or ALL
##'     \item ID: GO or KEGG IDs
##'     \item Description: description of functional category
##'     \item p.adjust: p-value adjusted for multiple hypothesis testing based 
##'     on method specified under pAdjustMethod argument
##'     \item qvalue: q value calculated with R’s qvalue function to control FDR
##'     \item itemID: IDs of items (genes for TSEA, drugs for DSEA) overlapping 
##'     among test and annotation sets.
##'     \item setSize: size of the functional category
##' } 
##' @slot organism organism information of the annotation system, 
##' only 'human' supported.
##' @slot ontology ontology type of the GO annotation system. If the annotation 
##' system is KEGG, it will be 'KEGG'
##' @slot drugs input drugs used for the enrichment test
##' @slot targets target information of the drugs in the defined drug-target 
##' annotation resource.
##' @exportClass feaResult
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