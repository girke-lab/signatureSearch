#' signatureSearch
#'
#' @name signatureSearch-package
#' @aliases signatureSearch-package signatureSearch gess fea 
#' @docType package
#' @useDynLib signatureSearch
#' @import Rcpp
#' @description 
#' Welcome to the signatureSearch package! This package integrates methods for
#' Gene Expression Signature Search (GESS) and Functional Enrichment Analyis 
#' (FEA). GESS methods use gene sets or Gene Expression Profiles (GEPs) 
#' of log2 fold change, z-scores or intensity value to search against large 
#' reference database containing tens to hundreds of thousands of GEPs or 
#' gene sets from compounds or genetic treatments. The large database was 
#' searched with very moderate memory requirements by setting chunks.
#' 
#' The FEA methods could be used to functionally annotate top ranking drugs in 
#' the GESS results in terms of GO categories or KEGG pathways.
#' @details 
#' The GESS methods include \code{CMAP}, \code{LINCS}, \code{gCMAP}, 
#' \code{Fisher} and \code{Cor}. For detailed
#' description, please see document of each method. Most methods 
#' can be easily parallelized for multiple query signatures.
#' 
#' The FEA can be done by two approaches, Target Set Enrichment Analysis (TSEA)
#' and Drug Set Enrichment Analysis (DSEA). The TSEA methods first get target
#' set of query drugs and do enrichment on the target protein/gene set. 
#' What should be noted is that different drugs may have the same targets. 
#' The target protein sets used for the TSEA contain often duplicated proteins. 
#' Standard enrichment methods would eliminate these duplications since they 
#' assume uniqueness in the test sets. Removing duplications in TSEA would be 
#' inappropriate since it would erase one of the most important pieces of 
#' information. To solve this problem, this project developed a weighting 
#' method for duplicated targets where the weighting is proportional to the 
#' frequency of the targets in the test set. The TSEA
#' includes \code{dup_hyperG}, \code{mGSEA} and \code{meanAbs} methods.
#' 
#' Instead of doing enrichment on the target set, the DSEA methods 
#' do enrichment directly on the drug set by mapping drugs to 
#' funcitonal categories via drug-target mapping information. Since drugs in 
#' the test sets are usually unique, it allows to use the classical enrichment
#' methods without any changes. The DSEA methods include \code{hyperG} and 
#' \code{GSEA}. For detailed description, please see document of each method.
#' @seealso 
#' Methods for GESS:
#'   \itemize{
#'     \item \code{\link{gess_cmap}}, \code{\link{gess_lincs}}, 
#'           \code{\link{gess_gcmap}} \code{\link{gess_fisher}}, 
#'           \code{\link{gess_cor}}
#'   }
#'    
#' Methods for FEA:
#'   \itemize{
#'      \item TSEA methods:    
#'          \code{\link{tsea_dup_hyperG}}, \code{\link{tsea_mGSEA}},
#'          \code{\link{tsea_mabs}}
#'      
#'      \item DSEA methods:
#'      \code{\link{dsea_hyperG}}, \code{\link{dsea_GSEA}}
#'   }
#' @author
#' \itemize{
#'   \item Yuzhu Duan (yduan004@ucr.edu)
#'   \item Thomas Girke (thomas.girke@ucr.edu)
#' } 
NULL



#' Drugs used in examples
#' 
#' A character vector containing top 10 drugs in the GESS result 
#' from gess_lincs method in the vignette
#'
#' @name drugs
#' @aliases drugs
#' @docType data
#' @keywords datasets
"drugs"

#' Target list used in examples
#' 
#' A named numeric vector with GENE SYMBOL as names. It is a subset of the 
#' first 1000 elements from 'targets' slot of 'mgsea_res' in the vignette. 
#' The scores represent the weights of targets/genes in the target set of the
#' selected top 10 drugs.
#'
#' @name targetList
#' @aliases targetList
#' @docType data
#' @keywords datasets
"targetList"