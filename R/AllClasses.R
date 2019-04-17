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
##' @slot refdb_name character(1), name of the reference database. Like "cmap",
##' "lincs" or other custom names.
##' @exportClass qSig
##' @keywords classes
setClass("qSig", slots = c(
  query = "ANY",
  gess_method = "character",
  refdb = "character",
  refdb_name = "character"
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
##'     \item pert: perturbation/drug names
##'     \item cell: cell types
##'     \item type: type of perturbation. In CMAP and LINCS databases, the 
##'     perturbation types are all treatment compound (trt_cp). Users can 
##'     build their custom signature database with other types of perturbation,
##'     e.g., gene knockdown or overexpression via \code{\link{build_custom_db}}
##'     function
##'     \item trend: up or down. up: the signature in the gess result is 
##'     positively connected with query signature; down: negatively connected.
##'     \item N_upset: number of genes in the up set of query signature
##'     \item N_downset: number of genes in the down set of query signature
##'     \item t_gn_sym: SYMBOL id of target genes/proteins of drugs.
##' } 
##' @slot query query signature
##' @slot gess_method method for GESS analysis
##' @slot refdb_name name of the reference database
##' @exportClass gessResult
##' @keywords classes
setClass("gessResult",
         slots = c(
           result = "data.frame",
           query = "ANY",
           gess_method = "character",
           refdb_name = "character"
         ))

## Constructor for "gessResult"
gessResult <- function(result, query, gess_method, refdb_name="UNKNOWN")
  new("gessResult", result=result, query=query, 
      gess_method=gess_method, refdb_name=refdb_name)


## Defining the validity method for "qSig"
# setValidity("qSig", function(object) {
#      TRUE
# })

##' Class "feaResult"
##' 
##' This class represents the result of functional enrichment analysis.
##'
##' The description of the tibble columns in the result slot:
##' \itemize{
##'     \item ID: GO term or KEGG pathway ID.
##'     \item Description: description of the functional categories
##'     \item GeneRatio: ratio of genes in the test set that are annotated
##'     at a specific GO node or KEGG pathway
##'     \item BgRatio: ratio of background genes that are annotated
##'     at a specific GO node or KEGG pathway
##'     \item pvalue: p value of the enrichment
##'     \item p.adjust: p value adjusted for multiple hypothesis testing using
#'     'p.adjust' function with defined method. 
##'     \item qvalue: q value
##'     \item geneID: genes/drugs overlapped between test set and annotation 
##'     sets.
##'     \item setSize: size of the functional category
##'     \item enrichmentScore: Enrichment Score (ES) from the GSEA algorithm 
##'     (Subramanian et al., 2005). It represents the degree to which a set S 
##'     is over-represented at the top or bottom of the scored ranked list L. 
##'     The score is calculated by walking down the list L, 
##'     increasing a running-sum statistic when we encounter a gene in S and 
##'     decreasing when it is not. The magnitude of the increment depends on 
##'     the gene scores (e.g., correlation of the gene with phenotype). 
##'     The ES is the maximum deviation from zero encountered in the random 
##'     walk; it corresponds to a weighted Kolmogorov-Smirnov-like statistic.
##'     \item NES: normalized enrichment score. The positive and negative
##'     enrichment scores are normalized separately by permutating the 
##'     gene labels of the gene list L 'nPerm' times and dividing the 
##'     enrichment score by mean of the permutaion ES with the same sign.
##'     \item pvalue: The nominal p-value of the ES is calculated using 
##'     permutation test. Specifically, the gene labels of the gene list L were
##'     permuted and the ES of the gene set was recomputed for the permutated 
##'     data, which generate a null distribution for the ES. The p-value of the 
##'     observed ES is then calculated relative to this null distribution.
##'     \item p.adjust: p values adjusted for multiple hypothesis testing 
##'     \item qvalues: q value calculated for FDR control 
##'     \item leadingEdge: genes in the gene set S (functional category) that 
##'     appear in the ranked list L at, or before, the point where the running 
##'     sum reaches its maximum deviation from zero. Can be interpreted as the 
##'     core of a gene set that accounts for the enrichment signal.
##'     \item ledge_rank: ranks of genes in 'leadingEdge' at the gene list L.
##'     \item mabs: Given a scored ranked gene list L, mabs(S) represents
##'     the mean absolute scores of the genes in set S. 
##'     \item Nmabs: In order to adjust for size variations in gene set S, 
##'     'nPerm' random permutations of L are performed to determine 
##'     permutation mabs. Subsequently, mabs(S) is normalized by subtracting 
##'     the median of the permutation mabs and then dividing by its standard 
##'     deviation yielding the normalized scores Nmabs(S).
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