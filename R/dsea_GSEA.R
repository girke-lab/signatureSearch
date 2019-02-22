##' Gene Set Enrichment Analysis of drug Ontology
##'
##'
##' @title dsea_GSEA
##' @param drugList ranked lists of all drugs in the GESS result. The scores of the corresponding GESS method can be used for ranking the drugs as it is 
##' required by the GSEA algorithm. The drugs with zero scores are excluded in the GESS results of CMAP, gCMAP and LINCS methods
##' @param type one of "GO" or "KEGG"
##' @param ont one of "BP", "MF", "CC" or "GO"
##' @param exponent weight of each step
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @export
##' @return A \code{feaResult} instance.
##' @author Yuzhu Duan
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
      # GO_DATA_drug <- get_GO_data_drug(OrgDb = "org.Hs.eg.db", ont, keytype="SYMBOL")
      # download GO_DATA_drug and save it locally to save time
      ext_path <- system.file("extdata", package="signatureSearch")
      godata_drug_path <- paste0(ext_path,"/GO_DATA_drug.rds")
      if(file.exists(godata_drug_path)){
        GO_DATA_drug <- readRDS(godata_drug_path)
      } else {
        download.file("http://biocluster.ucr.edu/~yduan004/fea/GO_DATA_drug.rds", godata_drug_path, quiet = TRUE)
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
      res@organism <- get_organism(OrgDb = "org.Hs.eg.db")
      res@ontology <- ont
      
      if (ont == "ALL") {
        res <- add_GO_Ontology(res, GO_DATA_drug)
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



