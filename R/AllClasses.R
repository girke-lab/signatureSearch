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

##' Class "gessResult"
##' 
##' The class stores the result of GESS analysis
##' @name gessResult-class
##' @aliases gessResult
##' @docType class
##' @slot result tibble from GESS analysis, represents a list of drugs 
##' in the reference database ranked by their signature similarity to the query.
##' 
##' The description of the common tibble columns from different GESS methods:
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
##' @slot gess_method method for GESS analysis
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

##' Class "feaResult"
##' 
##' This class represents the result of functional enrichment analysis.
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
##' @name feaResult-class
##' @aliases feaResult
##' @docType class
##' @slot result tibble representing enriched functional categories
##' @slot organism only "human" supported
##' @slot ontology biological ontology
##' @slot drugs Drug IDs
##' @slot targets Target IDs of drugs in DrugBank/LINCS/STITCH databases or 
##' target list with scores.
##' @exportClass feaResult
##' @author Yuzhu Duan
##' @references 
##' Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., 
##' Gillette, M. A., … Mesirov, J. P. (2005). Gene set enrichment analysis: a 
##' knowledge-based approach for interpreting genome-wide expression profiles. 
##' Proceedings of the National Academy of Sciences of the United States of
##' America, 102(43), 15545–15550. \url{https://doi.org/10.1073/pnas.0506580102}
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