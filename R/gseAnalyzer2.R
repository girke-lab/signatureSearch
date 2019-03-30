##' The modified Gene Set Enrichment Analysis (GSEA) of Gene Ontology supports
##' geneList with large portion of zeros and enrich GO terms at the top of the
##' geneList
##'
##' @title modified GSEA method for GO enrichment analysis
##' @param geneList numeric vector with gene SYMBOL as names, the genes are 
##' ranked decreasingly by their scores.
##' @param ont one of "BP", "MF", "CC" or "ALL"
##' @param OrgDb OrgDb, e.g., "org.Hs.eg.db".
##' @param keyType keytype of gene
##' @param exponent weight of each step
##' @param nproc If not equal to zero, sets BPPARAM to use nproc workers 
##' (default = 1)
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each gene set for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @return feaResult object
##' @examples 
##' \dontrun{
##' data(targetList)
##' library(org.Hs.eg.db)
##' gsego <- gseGO2(geneList=targetList, ont="MF", OrgDb=org.Hs.eg.db,
##'                 pvalueCutoff = 1)
##' head(gsego)
##' }
##' @export
gseGO2 <- function(geneList,
                  ont           = "BP",
                  OrgDb,
                  keyType       = "SYMBOL",
                  exponent      = 1,
                  nproc         = 1,
                  nPerm         = 1000,
                  minGSSize     = 2,
                  maxGSSize     = 500,
                  pvalueCutoff  = 0.05,
                  pAdjustMethod = "BH",
                  verbose       = TRUE) {

    ont %<>% toupper
    ont <- match.arg(ont, c("BP", "CC", "MF", "ALL"))
    #GO_DATA <- get_GO_data(OrgDb, ont, keytype="SYMBOL")
    # download GO_DATA.rds and save it to cache to save time
    fl <- download_data_file(url=
        "http://biocluster.ucr.edu/~yduan004/signatureSearch_data/GO_DATA.rds",
                             rname="GO_DATA")
    GO_DATA <- readRDS(fl)
    
    res <-  GSEA_internal2(geneList = geneList,
                          exponent = exponent,
                          nPerm = nPerm,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          verbose = verbose,
                          USER_DATA = GO_DATA,
                          nproc = nproc)

    if (is.null(res))
        return(res)
    # Add and select ontology in res
    res <- add_GO_Ontology(res, GO_DATA)
    tmp_df <- res@result
    colnames(tmp_df)[1] = "ont"
    res@result <- tmp_df
    if(ont != "ALL")
      res@result <- as_tibble(res[res$ont == ont, ])
    
    res@organism <- get_organism(OrgDb)
    return(res)
}

##' The modified GSEA algorithm supports geneList with large portion of zeros
##' and enrich KEGG pathways at the top of the geneList
##'
##' @title modified GSEA method for KEGG enrichment analysis
##' @param geneList order ranked geneList
##' @param organism supported organism listed in '
##' \url{http://www.genome.jp/kegg/catalog/org_list.html}
##' @param keyType one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
##' @param exponent weight of each step
##' @param nproc If not equal to zero, 
##' sets BPPARAM to use nproc workers (default = 1)
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @return feaResult object
##' @examples 
##' \dontrun{
##' # Gene Entrez id should be used for KEGG enrichment
##' data(geneList, package="DOSE")
##' geneList[100:length(geneList)]=0
##' gsekk <- gseKEGG2(geneList=geneList, pvalueCutoff = 1)
##' head(gsekk)
##' }
##' @export
gseKEGG2 <- function(geneList,
                    organism          = 'hsa',
                    keyType           = 'kegg',
                    exponent          = 1,
                    nproc             = 1,
                    nPerm             = 1000,
                    minGSSize         = 10,
                    maxGSSize         = 500,
                    pvalueCutoff      = 0.05,
                    pAdjustMethod     = "BH",
                    verbose           = TRUE) {
    species <- organismMapper(organism)
    KEGG_DATA <- prepare_KEGG(species, "KEGG", keyType)
    res <-  GSEA_internal2(geneList = geneList,
                          exponent = exponent,
                          nPerm = nPerm,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          verbose = verbose,
                          USER_DATA = KEGG_DATA,
                          nproc = nproc)

    if (is.null(res))
        return(res)

    res@organism <- species
    res@ontology <- "KEGG"

    return(res)
}



