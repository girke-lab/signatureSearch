setClassUnion("listOrMat", c("list", "matrix"))
setClassUnion("charOrNum", c("character", "numeric", "NULL"))

##' Class "qSig"
##' 
##' S4 object named \code{qSig} containing query signature information for Gene
##' Expression Signature (GES) searches. It contains slots for query signature,
##' GESS method and path to the GES reference database.
##' @name qSig-class
##' @docType class
##' @aliases qSig-class
##' @slot query If 'gess_method' is one of 'CMAP' or 'LINCS', 
##' this should be a list with two character vectors named \code{upset}
##' and \code{downset} for up- and down-regulated gene labels (here Entrez IDs),
##' repectively.
##' 
##' If 'gess_method' is 'gCMAP', 'Fisher' or 'Cor', a single column matrix with
##' gene expression values should be assigned. The corresponding gene labels 
##' are stored in the row name slot of the matrix. The expected type of gene 
##' expression values is explained in the help files of the corresponding GESS 
##' methods.
##' @slot gess_method one of 'CMAP', 'LINCS', 'gCMAP', 'Fisher' or 'Cor'
##' @slot refdb character(1), can be "cmap", "cmap_expr", "lincs", or 
##' "lincs_expr" when using existing CMAP/LINCS databases. 
##' 
##' If users want to use a custom signature database, it should be the file path
##' to the HDF5 file generated with the \code{\link{build_custom_db}} function.
##' Alternatively, source files of the CMAP/LINCS databases can be used as
##' explained in the vignette of the
##' \code{\link[signatureSearchData]{signatureSearchData}} package.
##' @exportClass qSig
##' @keywords classes
setClass("qSig", slots = c(
  query="listOrMat",
  gess_method = "character",
  refdb = "character"
))

##' gessResult object
##' 
##' The \code{gessResult} object organizes Gene Expression Signature Search 
##' (GESS) results. This includes the main tabular result of a GESS, its query 
##' signature, the name of the chosen GESS method and the path to the reference
##' database.
##' @name gessResult-class
##' @aliases gessResult-class
##' @docType class
##' @slot result tibble object containing the search results for each
##' perturbagen (e.g. drugs) in the reference database ranked by their
##' signature similarity to the query. The result table can be extracted via
##' the \code{\link{result}} accessor function.
##' 
##' Descriptions of the columns common among the tabular results of the 
##' individual GESS methods are given below. Note, the columns specific to each 
##' GESS method are described in their help files.
##' \itemize{
##'     \item pert: character, name of perturbagen (e.g. drug) in the reference 
##'     database
##'     \item cell: character, acronym of cell type
##'     \item type: character, perturbation type. In the CMAP/LINCS 
##'     databases provided by \code{signatureSearchData}, the perturbation types
##'     are currently treatments with drug-like compounds (trt_cp). If required,
##'     users can build custom signature database with other types of
##'     perturbagens (e.g., gene knockdown or over-expression events) with the 
##'     provided \code{\link{build_custom_db}} function.
##'     \item trend: character, up or down when the reference signature is 
##'     positively or negatively connected with the query signature, 
##'     respectively.
##'     \item N_upset: integer, number of genes in the query up set
##'     \item N_downset: integer, number of genes in the query down set
##'     \item t_gn_sym: character, symbol of the gene encoding the
##'     corresponding drug target protein
##' } 
##' @slot query query signature
##' @slot gess_method name of the GESS method 
##' @slot refdb path to the reference database
##' @exportClass gessResult
##' @keywords classes
setClass("gessResult",
         slots = c(
           result = "data.frame",
           query = "listOrMat",
           gess_method = "character",
           refdb = "character"
         ))

#' This is a helper function to construct a \code{gessResult} object. For 
#' detail description, please consult the help file of the 
#' \code{\link{gessResult-class}}.
#' @title Constructor for \code{\link{gessResult-class}}
#' @param result tibble object containing the GESS results
#' @param query list or a matrix, query signature
#' @param gess_method character(1), name of the GESS method
#' @param refdb character(1), path to the reference database
#' @return \code{gessResult} object
#' @examples 
#' gr <- gessResult(result=dplyr::tibble(pert=letters[seq_len(10)], 
#'                                       val=seq_len(10)), 
#'                  query=list(up=c("g1","g2"), down=c("g3","g4")),
#'                  gess_method="LINCS", refdb="path/to/lincs/db")
#' @export 
gessResult <- function(result, query, gess_method, refdb)
  new("gessResult", result=result, query=query, 
      gess_method=gess_method, refdb=refdb)


## Defining the validity method for "qSig"
# setValidity("qSig", function(object) {
#      TRUE
# })

##' feaResult object
##' 
##' The \code{feaResult} object stores Functional Enrichment Analysis (FEA) 
##' results generated by the corresponding Target and Drug Set Enrichment 
##' methods (here TSEA and DSEA) defined by \code{signatureSearch}. This 
##' includes slots for the FEA results in tabular format, the organism 
##' information, and the type of functional annotation used (e.g. GO or KEGG). 
##' It also includes the drug information used for the FEA, as well as the 
##' corresponding target protein information.
##' @name feaResult-class
##' @aliases feaResult-class
##' @docType class
##' @slot result tibble object, this tabular result contains the
##' enriched functional categories (e.g. GO terms or KEGG pathways) ranked by 
##' the corresponding enrichment statistic. The result table can be extracted 
##' via the \code{\link{result}} accessor function.
##' 
##' Description of the columns that are shared among the result tables 
##' generated by the different FEA methods:
##' \itemize{
##'     \item ont: in case of GO, one of BP, MF, CC, or ALL
##'     \item ID: GO or KEGG IDs
##'     \item Description: description of functional category
##'     \item p.adjust: p-value adjusted for multiple hypothesis testing based 
##'     on method specified under pAdjustMethod argument
##'     \item qvalue: q value calculated with R's qvalue function to control FDR
##'     \item itemID: IDs of items (genes for TSEA, drugs for DSEA) overlapping 
##'     among test and annotation sets.
##'     \item setSize: size of the functional category
##' } 
##' @slot organism organism information of the annotation system. 
##' Currently, limited to 'human', since drug-target annotations are too sparse
##' for other organisms.
##' @slot ontology ontology type of the GO annotation system. If the 
##' annotation system is KEGG, it will be 'KEGG'
##' @slot drugs input drug names used for the enrichment test
##' @slot targets target information for the query drugs obtained from the 
##' chosen drug-target annotation source.
##' @exportClass feaResult
##' @keywords classes
setClass("feaResult",
         slots = c(
           result         = "data.frame",
           organism       = "character",
           ontology       = "character",
           drugs          = "character",
           targets        = "charOrNum"
         )
)

#' This is a helper function to construct a \code{feaResult} object. For 
#' detail description, please consult the help file of the 
#' \code{\link{feaResult-class}}.
#' @title Constructor for \code{\link{feaResult-class}}
#' @param result tibble object containing the FEA results
#' @param organism character(1), organism information of the annotation system
#' @param ontology character(1), ontology type of the GO annotation system. 
#' If the annotation system is KEGG, it will be 'KEGG'
#' @param drugs character vector, input drug names used for the enrichment test
#' @param targets character vector, gene labels of the gene/protein targets 
#' for the drugs 
#' @return \code{feaResult} object
#' @examples
#' fr <- feaResult(result=dplyr::tibble(id=letters[seq_len(10)], 
#'                                      val=seq_len(10)),
#'                 organism="human", ontology="MF", drugs=c("d1", "d2"), 
#'                 targets=c("t1","t2"))
#' @export 
feaResult <- function(result, organism="UNKNOWN", ontology="UNKNOWN", 
                      drugs="UNKNOWN", targets="UNKNOWN")
  new("feaResult", result=result, organism=organism, ontology=ontology,
      drugs=drugs, targets=targets)
