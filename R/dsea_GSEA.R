##' The original GSEA method is used to do enrichment analysis on a scored
##' ranked drug list after mapping drugs to functional categories via 
##' drug-target links in DrugBank, CLUE and STITCH databases
##' 
##' The ranked drug list consists of all drugs in the GESS result. 
##' The similarity scores of the corresponding GESS method can be used for 
##' ranking the drugs as it is required by the GSEA algorithm. The drugs 
##' with zero scores are excluded. 
##'
##' @title GSEA method for DSEA
##' @param drugList scored ranked list of all drugs in the GESS result. 
##' The similarity scores of the corresponding GESS method can be used for 
##' ranking the drugs as it is required by the GSEA algorithm. 
##' The drugs with zero scores are excluded.
##' @param type one of "GO" or "KEGG"
##' @param ont one of "BP", "MF", "CC" or "GO"
##' @param exponent weight of each step
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of drug sets annotated by ontology term 
##' after drug to functional category mappings.
##' @param maxGSSize maximal size of drug sets annotated for testing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod p value adjustment method,
##' one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @return \code{\link{feaResult}} object, 
##' represents enriched functional categories.
##' @seealso \code{\link{feaResult}}, \code{\link{fea}},
##'          \code{\link[signatureSearch_data]{dtlink_db_clue_sti.db}}
##' @references 
##' GSEA algorithm: Subramanian et at., 2005,
##' \url{https://www.pnas.org/content/102/43/15545}
##' @examples 
##' db_dir <- system.file("extdata", "sample_db", package = "signatureSearch")
##' sample_db <- loadHDF5SummarizedExperiment(db_dir)
##' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
##' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
##' query = as.numeric(query_mat); names(query) = rownames(query_mat)
##' upset <- head(names(query[order(-query)]), 150)
##' downset <- tail(names(query[order(-query)]), 150)
##' qsig_lincs <- qSig(qsig = list(upset=upset, downset=downset), 
##'                    gess_method = "LINCS", refdb = sample_db)
##' lincs <- gess_lincs(qsig_lincs, sortby="NCS")
##' dl <- abs(result(lincs)$NCS); names(dl) <- result(lincs)$pert
##' dl <- dl[dl>0]
##' dl <- dl[!duplicated(names(dl))]
##' # GO annotation system
##' gsea_res <- dsea_GSEA(drugList=dl, type="GO", ont="MF", exponent=1, 
##'                       nPerm=1000, pvalueCutoff=0.2, minGSSize=5)
##'                       result(gsea_res)
##' # KEGG annotation system
##' gsea_k_res <- dsea_GSEA(drugList=dl, type="KEGG", exponent=1, nPerm=1000, 
##'                         pvalueCutoff=0.5, minGSSize=5)
##' result(gsea_k_res)
##' @export
dsea_GSEA <- function(drugList, 
                      type = "GO", 
                      ont           = "BP", 
                      exponent      = 1,
                      nPerm         = 1000,
                      minGSSize     = 10,
                      maxGSSize     = 500,
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH") {
  
  ont %<>% toupper
  ont <- match.arg(ont, c("BP", "CC", "MF", "ALL"))
  names(drugList) <- tolower(names(drugList))
  # check whether there are duplicated names in drugList
  if(sum(duplicated(names(drugList))) > 0)
    stop("The names of the query drug list for dsea_GSEA need to be unique!")
  
  if(type=="GO"){
    # GO_DATA_drug <- get_GO_data_drug(OrgDb = "org.Hs.eg.db", 
    #                                  ont, keytype="SYMBOL")
    # download GO_DATA_drug and save it locally to save time
    ext_path <- system.file("extdata", package="signatureSearch")
    godata_drug_path <- paste0(ext_path,"/GO_DATA_drug.rds")
    if(file.exists(godata_drug_path)){
      GO_DATA_drug <- readRDS(godata_drug_path)
    } else {
      download.file("http://biocluster.ucr.edu/~yduan004/fea/GO_DATA_drug.rds", 
                    godata_drug_path, quiet = TRUE)
      GO_DATA_drug <- readRDS(godata_drug_path)
    }
    
    res <-  GSEA_internal(geneList = drugList,
                          exponent = exponent,
                          nPerm = nPerm,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          USER_DATA = GO_DATA_drug,
                          seed = FALSE)
    
    if (is.null(res))
      return(res)
    res@organism <- clusterProfiler:::get_organism(OrgDb = "org.Hs.eg.db")
    res@ontology <- ont
    
    if (ont == "ALL") {
      res <- clusterProfiler:::add_GO_Ontology(res, GO_DATA_drug)
    } 
    return(res)
  }
  
  if(type=="KEGG"){
    species <- organismMapper("hsa")
    KEGG_DATA_drug <- prepare_KEGG_drug(species, "KEGG", keyType="kegg")
    
    res <-  GSEA_internal(geneList = drugList,
                          exponent = exponent,
                          nPerm = nPerm,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          USER_DATA = KEGG_DATA_drug,
                          seed = FALSE)
    
    if (is.null(res))
      return(res)
    res@organism <- species
    res@ontology <- "KEGG"
    return(res)
  }
}



