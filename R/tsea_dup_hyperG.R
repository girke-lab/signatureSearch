#' dup_hyperG method for TSEA
#' 
#' This method support target set with duplications by adjusting the frequency 
#' of the duplicated proteins in the target set instead of taking the unique.
#' 
#' The classical hypergeometric test assumes uniqueness in its gene/protein 
#' test sets. To maintain the duplication information in the test sets used for 
#' TSEA, the duplication information in the test set is maintained by adjusting
#' their frequency.
#' @param drugs query drug set used to do TSEA.
#' Can be top ranking drugs in the GESS result. 
#' @param universe background genes/targets. If set as `Default`, 
#' it represents all the annotated genes/targets in
#' the corresponding annotation system (e.g. GO or KEGG). 
#' If "type" is "GO", universe should be gene SYMBOL ids. If "type" is "KEGG",
#' universe should be gene entrez ids.
#' @param type one of `GO` or `KEGG`
#' @param ont if type is `GO`, set ontology. One of `BP`,`MF`,`CC` or `ALL`
#' @param pAdjustMethod p value adjustment method, 
#' one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param pvalueCutoff p value Cutoff
#' @param qvalueCutoff qvalue Cutoff
#' @param minGSSize minimum size of each gene set in the annotation system
#' @param maxGSSize maximum size of each gene set in the annotation system
#' @return \code{\link{feaResult}} object, 
#' represents enriched functional categories.
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
tsea_dup_hyperG <- function(drugs, universe = "Default", 
                            type="GO", ont="MF", 
                            pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.05, 
                            minGSSize = 5, maxGSSize = 500){
  # message("The query drugs [", length(drugs),"] are: \n", 
  #         paste0(drugs[seq_len(min(length(drugs), 10))], sep ="  "), "...")
  drugs <- unique(tolower(drugs))
  targets <- get_targets(drugs, database = "all")
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