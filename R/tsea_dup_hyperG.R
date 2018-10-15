#' drug targets GO/KEGG enrichment analysis by using `dup_hyperG` test or modified `gsea` method
#' @title tsea_dup_hyperG
#' @param drugs query drug set used to do target set enrichment analysis
#' @param universe background genes/targets. If set as `Default`, it represents all the annotated genes/targets in
#' the corresponding annotation system. If "type" is "GO", universe should be gene SYMBOL ids. If "type" is "KEGG",
#' universe should be gene entrez ids.
#' @param type can be `GO` or `KEGG`
#' @param ont if type is `GO`, set ontology, can be `BP`,`MF`,`CC` or `ALL`
#' @param pAdjustMethod p value adjustment method for p values in target set enrichment result
#' @param pvalueCutoff p value Cutoff
#' @param qvalueCutoff qvalue Cutoff
#' @param minGSSize minimum size of each gene set in annotation system
#' @param maxGSSize maximum size of each gene set in annotation system
#' @return feaResult object
#' @export
#' @author Yuzhu Duan (yduan004@ucr.edu)
#' 
tsea_dup_hyperG <- function(drugs, universe = "Default", 
                      type="GO", ont="MF", 
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
                      minGSSize = 2, maxGSSize = 500){
  drugs <- unique(tolower(drugs))
  message("The query drugs [", length(drugs),"] are: \n", paste0(drugs, collapse =" \ "))
  targets <- get_targets(drugs, database = "all")
  gnset <- na.omit(unlist(lapply(targets$t_gn_sym, function(i) unlist(strsplit(as.character(i), split = ";")))))
  
  # tmpfile <- tempfile(fileext=".txt")
  # download.file("http://biocluster.ucr.edu/~yduan004/db_lincs_sti_univ_gnset.txt", tmpfile, quiet = TRUE)
  # universe <- readLines(tmpfile)
  # unlink(tmpfile)
  
  if(type=="GO"){
    if(universe=="Default"){
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
    }
    ego <- enrichGO(gene = gnset, universe = universe, OrgDb = org.Hs.eg.db, keytype = 'SYMBOL', 
                                      ont = ont, pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, 
                                      qvalueCutoff = qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize)
    ego@drugs=drugs
    return(ego)
  }
  if(type=="KEGG"){
    if(universe=="Default"){
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
    }
    #univ_entrez <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = universe, keytype = "SYMBOL", columns = "ENTREZID")$ENTREZID)
    gnset_entrez <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = gnset, keytype = "SYMBOL", columns = "ENTREZID")$ENTREZID)
    kk <- enrichKEGG(gene=gnset_entrez, organism = "hsa", keyType = "kegg", pvalueCutoff = pvalueCutoff, 
                                       pAdjustMethod = pAdjustMethod, universe=universe, minGSSize = minGSSize, 
                                       maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
    kk@drugs <- drugs
    return(kk)
  }
}