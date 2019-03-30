##' MeanAbs Analysis of Gene Ontology
##'
##' @title mabsGO
##' @param geneList scored ranked geneList
##' @param ont one of "BP", "MF", "CC" or "ALL"
##' @param OrgDb OrgDb
##' @param keyType keytype of gene
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @return \code{\link{feaResult}} object
##' @author Yuzhu Duan
##' @examples 
##' \dontrun{
##' data(targetList)
##' library(org.Hs.eg.db)
##' mg <- mabsGO(geneList=targetList, ont="MF", OrgDb=org.Hs.eg.db,
##'              pvalueCutoff = 1)
##' head(mg)
##' }
##' @export
mabsGO <- function(geneList,
                  ont           = "BP",
                  OrgDb,
                  keyType       = "SYMBOL",
                  nPerm         = 1000,
                  minGSSize     = 5,
                  maxGSSize     = 500,
                  pvalueCutoff  = 0.05,
                  pAdjustMethod = "BH") {

    ont %<>% toupper
    ont <- match.arg(ont, c("BP", "CC", "MF", "ALL"))

    #GO_DATA <- get_GO_data(OrgDb, ont, keytype="SYMBOL")
    # download GO_DATA.rds and save it to cache to save time
    fl <- download_data_file(url=
        "http://biocluster.ucr.edu/~yduan004/signatureSearch_data/GO_DATA.rds",
                             rname="GO_DATA")
    GO_DATA <- readRDS(fl)
    
    res <-  mabs_internal(geneList = geneList,
                          nPerm = nPerm,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          USER_DATA = GO_DATA)

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
    res@ontology <- ont
    return(res)
}

##' meanAbs analysis of KEGG
##'
##' @title mabsKEGG
##' @param geneList scored ranked geneList
##' @param organism supported organism listed in 
##' \url{http://www.genome.jp/kegg/catalog/org_list.html}
##' @param keyType one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @return \code{\link{feaResult}} object
##' @examples 
##' \dontrun{
##' # Gene Entrez id should be used for KEGG enrichment
##' data(geneList, package="DOSE")
##' geneList[100:length(geneList)]=0
##' mk <- mabsKEGG(geneList=geneList, pvalueCutoff = 1)
##' head(mk)
##' }
##' @export
mabsKEGG <- function(geneList,
                    organism          = 'hsa',
                    keyType           = 'kegg',
                    nPerm             = 1000,
                    minGSSize         = 5,
                    maxGSSize         = 500,
                    pvalueCutoff      = 0.05,
                    pAdjustMethod     = "BH") {

    species <- organismMapper(organism)
    KEGG_DATA <- prepare_KEGG(species, "KEGG", keyType)

    res <-  mabs_internal(geneList = geneList,
                          nPerm = nPerm,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          USER_DATA = KEGG_DATA)

    if (is.null(res))
        return(res)

    res@organism <- species
    res@ontology <- "KEGG"
    return(res)
}

