#' dup_hyperG Enrichment Method
#' 
#' This duplication adjusted hypergeometric test supports target gene/protein 
#' set with duplications by adjusting 
#' the frequency of the duplicated proteins in the target set instead of 
#' taking the unique value.
#' @details 
#' The classical hypergeometric test assumes uniqueness in its gene/protein 
#' test sets. To maintain the duplication information in the test sets used 
#' for TSEA, the values of total number of genes in the test set and the number
#' of genes in the test set annotated at a functional category are adjusted by 
#' maintating the frequency of the target proteins in the test set instead of 
#' taking the unique value.
#' @section Column description:
#' Description of the columns in the result table specific to the 
#' hypergeometric test: 
#' \itemize{
#'     \item GeneRatio: ratio of genes in the test set that are annotated at a 
#'     specific GO node or KEGG pathway
#'     \item BgRatio: ratio of background genes that are annotated
#'     at a specific GO node or KEGG pathway
#'     \item pvalue: raw p-value of enrichment test
#' }
#' Description of the other columns are available at the 'result' slot of the
#' \code{\link{feaResult}} object.
#' 
#' @param drugs character vector, query drug set used for functional enrichment.
#' Can be top ranking drugs in the GESS result. 
#' @param universe character vector, labels of background genes/proteins. 
#' If set as 'Default', it uses all the annotated genes/proteins in
#' the corresponding annotation system (e.g. GO or KEGG). 
#' If 'type' is 'GO', it should be a vector of gene SYMBOL IDs. 
#' If 'type' is 'KEGG', it should be gene Entrez IDs.
#' @param type one of `GO` or `KEGG`
#' @param ont character(1). If type is `GO`, set ontology as one of `BP`,`MF`,
#' `CC` or `ALL`. If type is 'KEGG', it is ignored.
#' @param pAdjustMethod p-value adjustment method, 
#' one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'
#' @param pvalueCutoff double, p-value cutoff
#' @param qvalueCutoff double, qvalue cutoff
#' @param minGSSize integer, minimum size of each gene set in annotation system
#' @param maxGSSize integer, maximum size of each gene set in annotation system
#' @param dt_anno drug-target annotation resource. one of 'DrugBank', 'CLUE', 
#' 'STITCH' or 'all'. If 'dt_anno' is 'all', the targets from DrugBank, CLUE 
#' and STITCH databases will be combined. It is recommended to set the 'dt_anno'
#' as 'all' since it will get the most complete target set of as many drugs
#' as possible. Users could also choose individual annotation resource if 
#' wanted, but should be aware that if the chosen drug-target annotation
#' resource contains limited drugs (such as CLUE), many query drugs will not
#' get targets and the target set may not be complete, which could affect
#' the enrichment result.  
#' @return \code{\link{feaResult}} object, the result table contains the
#' enriched functional categories (e.g. GO terms or KEGG pathways) ranked by 
#' the corresponding enrichment statistic.
#' @seealso \code{\link{feaResult}}, \code{\link{fea}}
#' @examples 
#' data(drugs)
#' \dontrun{
#' ## GO annotation system
#' dup_hyperG_res <- tsea_dup_hyperG(drugs = drugs, universe = "Default", 
#'                                   type = "GO", ont="MF", pvalueCutoff=0.05,
#'                                   pAdjustMethod="BH", qvalueCutoff = 0.1, 
#'                                   minGSSize = 10, maxGSSize = 500)
#' result(dup_hyperG_res)
#' }
#' ## KEGG annotation system
#' dup_hyperG_k_res <- tsea_dup_hyperG(drugs = drugs, universe = "Default", 
#'                                     type = "KEGG", pvalueCutoff=0.1, 
#'                                     pAdjustMethod="BH", qvalueCutoff = 0.2, 
#'                                     minGSSize = 10, maxGSSize = 500)
#' result(dup_hyperG_k_res)
#' @export
tsea_dup_hyperG <- function(drugs, universe="Default", 
                            type="GO", ont="MF", 
                            pAdjustMethod="BH", pvalueCutoff=0.05, 
                            qvalueCutoff=0.05, 
                            minGSSize=5, maxGSSize=500,
                            dt_anno="all"){
  # message("The query drugs [", length(drugs),"] are: \n", 
  #         paste0(drugs[seq_len(min(length(drugs), 10))], sep ="  "), "...")
  drugs <- unique(tolower(drugs))
  targets <- get_targets(drugs, database = dt_anno)
  gnset <- na.omit(unlist(lapply(targets$t_gn_sym, function(i) 
    unlist(strsplit(as.character(i), split = "; ")))))
  
  if(type=="GO"){
    if(universe=="Default"){
      # Get universe genes in GO annotation system
      ext_path <- system.file("extdata", package="signatureSearch")
      univ_tar_path <- paste0(ext_path,"/univ_genes_in_GO.txt")
      universe <- readLines(univ_tar_path)
    }
    ego <- enrichGO2(gene = gnset, universe = universe, OrgDb = org.Hs.eg.db, 
                    keytype = 'SYMBOL', ont = ont, 
                    pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, 
                    qvalueCutoff = qvalueCutoff, 
                    minGSSize = minGSSize, maxGSSize = maxGSSize)
    ego@drugs=drugs
    return(ego)
  }
  if(type=="KEGG"){
    if(universe=="Default"){
      # Get universe genes in KEGG annotation system (entrez ids)
      ext_path <- system.file("extdata", package="signatureSearch")
      univ_tar_path <- paste0(ext_path,"/univ_genes_in_KEGG.txt")
      universe <- readLines(univ_tar_path)
    }
    gnset_entrez <- suppressMessages(
      AnnotationDbi::select(org.Hs.eg.db, keys = gnset, keytype = "SYMBOL", 
                            columns = "ENTREZID")$ENTREZID)
    kk <- enrichKEGG2(gene=gnset_entrez, organism = "hsa", keyType = "kegg", 
                     pvalueCutoff = pvalueCutoff, 
                     pAdjustMethod = pAdjustMethod, 
                     universe=universe, minGSSize = minGSSize, 
                     maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
    kk@drugs <- drugs
    return(kk)
  }
}