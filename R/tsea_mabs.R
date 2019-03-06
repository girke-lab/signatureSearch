#' tsea_mabs
#' 
#' drug targets GO/KEGG enrichment analysis by using `mabs` method
#' @param drugs query drug set used to do TSEA, Can be top ranking drugs in GESS result. 
#' @param type can be `GO` or `KEGG`
#' @param ont if type is `GO`, set ontology, can be `BP`,`MF`,`CC` or `ALL`
#' @param nPerm permutation numbers used to calculate p value
#' @param pAdjustMethod p value adjustment method for p values in the enrichment result
#' @param pvalueCutoff p value Cutoff
#' @param minGSSize minimum size of each gene set in annotation system
#' @param maxGSSize maximum size of each gene set in annotation system
#' @return A \code{feaResult} instance
#' @export
#' @author Yuzhu Duan (yduan004@ucr.edu)
#' 
tsea_mabs <- function(drugs, 
                      type="GO", ont="MF", 
                      nPerm=1000,  
                      pAdjustMethod = "BH", pvalueCutoff = 0.05,
                      minGSSize = 5, maxGSSize = 500){
  drugs <- unique(tolower(drugs))
  targets <- get_targets(drugs, database = "all")
  gnset <- na.omit(unlist(lapply(targets$t_gn_sym, function(i) 
    unlist(strsplit(as.character(i), split = "; ")))))
  # give scores to gnset
  tar_tab <- table(gnset)
  tar_dup <- as.numeric(tar_tab); names(tar_dup) <- names(tar_tab)
  tar_weight <- sort(tar_dup/sum(tar_dup), decreasing = TRUE)
  
  if(type=="GO"){
    # Get universe genes in GO annotation system
    ext_path <- system.file("extdata", package="signatureSearch")
    univ_tar_path <- paste0(ext_path,"/univ_genes_in_GO.txt")
    universe <- readLines(univ_tar_path)
    tar_diff <- setdiff(universe, gnset)
    tar_diff_weight <- rep(0, length(tar_diff))
    names(tar_diff_weight)=tar_diff
    tar_total_weight = c(tar_weight, tar_diff_weight)
    mabsgo <- mabsGO(geneList = tar_total_weight, OrgDb = org.Hs.eg.db, 
                     ont = ont, keyType = "SYMBOL", nPerm = nPerm, 
                     minGSSize = minGSSize, maxGSSize = maxGSSize, 
                     pvalueCutoff = pvalueCutoff, pAdjustMethod=pAdjustMethod)
    mabsgo@drugs = drugs
    return(mabsgo)
  }
  
  if(type=="KEGG"){
    # Get universe genes in KEGG annotation system
    ext_path <- system.file("extdata", package="signatureSearch")
    univ_tar_path <- paste0(ext_path,"/univ_genes_in_KEGG.txt")
    universe <- readLines(univ_tar_path)
    # convert gnset symbol to entrez
    gnset_entrez <- suppressMessages(
      AnnotationDbi::select(org.Hs.eg.db, keys = gnset, 
                            keytype = "SYMBOL", columns = "ENTREZID"))
    gnset_entrez2 <- as.character(na.omit(gnset_entrez$ENTREZID))
    # give scores to gnset_entrez2
    tar_tab <- table(gnset_entrez2)
    tar_dup <- as.numeric(tar_tab); names(tar_dup) <- names(tar_tab)
    tar_weight <- sort(tar_dup/sum(tar_dup), decreasing = TRUE)
    
    tar_diff <- setdiff(universe, gnset_entrez2)
    tar_diff_weight <- rep(0, length(tar_diff)); names(tar_diff_weight)=tar_diff
    tar_total_weight = c(tar_weight, tar_diff_weight)

    mabskk <- mabsKEGG(geneList=tar_total_weight, organism='hsa', 
                       keyType='kegg', nPerm = nPerm, 
                       minGSSize = minGSSize, maxGSSize=maxGSSize, 
                       pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod)
    if(is.null(mabskk))
      return(NULL)
    mabskk@drugs = drugs
    return(mabskk)
  }
}