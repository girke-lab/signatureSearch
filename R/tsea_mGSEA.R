
#' drug targets GO/KEGG enrichment analysis by using `mGSEA` method
#' @title tsea_mGSEA
#' @param drugs query drug set used to do target set enrichment analysis (TSEA), Can be top ranking drugs in GESS result. 
#' @param type can be `GO` or `KEGG`
#' @param ont if type is `GO`, set ontology, can be `BP`,`MF`,`CC` or `ALL`
#' @param nPerm permutation numbers used to calculate p value
#' @param exponent weight of each step
#' @param pAdjustMethod p value adjustment method for p values in TSEA result
#' @param pvalueCutoff p value Cutoff
#' @param minGSSize minimum size of each gene set in annotation system
#' @param maxGSSize maximum size of each gene set in annotation system
#' @return A \code{feaResult} instance
#' @export
#' @author Yuzhu Duan (yduan004@ucr.edu)
#' 
tsea_mGSEA <- function(drugs, 
                      type="GO", ont="MF", 
                      nPerm=1000, exponent=1, 
                      pAdjustMethod = "BH", pvalueCutoff = 0.05,
                      minGSSize = 2, maxGSSize = 500){
  drugs <- unique(tolower(drugs))
  #message("The query drugs [", length(drugs),"] are: \n", paste0(drugs, collapse =" \ "))
  targets <- get_targets(drugs, database = "all")
  gnset <- na.omit(unlist(lapply(targets$t_gn_sym, function(i) unlist(strsplit(as.character(i), split = ";")))))
  # give scores to gnset
  tar_tab <- table(gnset)
  tar_dup <- as.numeric(tar_tab); names(tar_dup) <- names(tar_tab)
  tar_weight <- sort(tar_dup/sum(tar_dup), decreasing = T)
  
  if(type=="GO"){
    # Get universe genes in GO annotation system
    ext_path <- system.file("extdata", package="signatureSearch")
    univ_tar_path <- paste0(ext_path,"/univ_genes_in_GO.txt")
    if(file.exists(univ_tar_path)){
      universe <- readLines(univ_tar_path)
    } else {
      tryCatch(download.file("http://biocluster.ucr.edu/~yduan004/fea/univ_genes_in_GO.txt", univ_tar_path, quiet = TRUE), 
               error = function(e){
                 stop("Error happens when downloading fea/univ_genes_in_GO.txt")
               }, warning = function(w) {file.remove(univ_tar_path)
                 stop("Error happens when downloading fea/univ_genes_in_GO.txt")})
      universe <- readLines(univ_tar_path)
    }
    tar_diff <- setdiff(universe, gnset)
    tar_diff_weight <- rep(0, length(tar_diff)); names(tar_diff_weight)=tar_diff
    tar_total_weight = c(tar_weight, tar_diff_weight)
    gsego <- gseGO2(geneList = tar_total_weight, OrgDb = org.Hs.eg.db, ont = ont, keyType = "SYMBOL", 
                    nPerm = nPerm, minGSSize = minGSSize, maxGSSize = maxGSSize, pvalueCutoff = pvalueCutoff, 
                    exponent = exponent, nproc=1, verbose=TRUE, pAdjustMethod=pAdjustMethod)
    gsego@drugs = drugs
    return(gsego)
  }
  
  if(type=="KEGG"){
    # Get universe genes in KEGG annotation system
    ext_path <- system.file("extdata", package="signatureSearch")
    univ_tar_path <- paste0(ext_path,"/univ_genes_in_KEGG.txt")
    if(file.exists(univ_tar_path)){
      universe <- readLines(univ_tar_path)
    } else {
      tryCatch(download.file("http://biocluster.ucr.edu/~yduan004/fea/univ_genes_in_KEGG.txt", univ_tar_path, quiet = TRUE), 
               error = function(e){
                 stop("Error happens when downloading fea/univ_genes_in_KEGG.txt")
               }, warning = function(w) {file.remove(univ_tar_path)
                 stop("Error happens when downloading fea/univ_genes_in_KEGG.txt")})
      universe <- readLines(univ_tar_path)
    }
    
    tar_diff <- setdiff(universe, gnset)
    tar_diff_weight <- rep(0, length(tar_diff)); names(tar_diff_weight)=tar_diff
    tar_total_weight = c(tar_weight, tar_diff_weight)
    
    sym_entrez <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = names(tar_total_weight), keytype = "SYMBOL", columns = "ENTREZID"))
    ttw_entr <- NULL
    for(i in seq_along(tar_total_weight)){
      item = tar_total_weight[i]
      if(length(sym_entrez[sym_entrez$SYMBOL==names(item),"ENTREZID"])>0){
        names(item) <- sym_entrez[sym_entrez$SYMBOL==names(item),"ENTREZID"][1]
        ttw_entr = c(ttw_entr, item)
      }
    }
    gsekk <- gseKEGG2(geneList=ttw_entr, organism='hsa', keyType='kegg', nPerm = nPerm, 
                                       minGSSize = minGSSize, maxGSSize=maxGSSize, pvalueCutoff=pvalueCutoff, 
                                       pAdjustMethod = pAdjustMethod, verbose = TRUE)
    gsekk@drugs = drugs
    return(gsekk)
  }
}
