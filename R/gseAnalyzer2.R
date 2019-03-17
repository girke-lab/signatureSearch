##' Modified Gene Set Enrichment Analysis (GSEA) of Gene Ontology.
##' The GSEA method is modified to accept geneList with large portion of zeros
##'
##' @title gseGO2
##' @param geneList scored ranked geneList
##' @param ont one of "BP", "MF", "CC" or "ALL"
##' @param OrgDb OrgDb
##' @param keyType keytype of gene
##' @param exponent weight of each step
##' @param nproc If not equal to zero, sets BPPARAM to use nproc workers (default = 1)
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @export
##' @return feaResult object
gseGO2 <- function(geneList,
                  ont           = "BP",
                  OrgDb,
                  keyType       = "ENTREZID",
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
    # download GO_DATA and save it locally to save time
    ext_path <- system.file("extdata", package="signatureSearch")
    godata_path <- paste0(ext_path,"/GO_DATA.rds")
    if(file.exists(godata_path)){
      GO_DATA <- readRDS(godata_path)
    } else {
      download.file("http://biocluster.ucr.edu/~yduan004/fea/GO_DATA.rds", godata_path, quiet = TRUE)
      GO_DATA <- readRDS(godata_path)
    }
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

    res@organism <- get_organism(OrgDb)

    if (ont == "ALL") {
        res <- add_GO_Ontology(res, GO_DATA)
    }
    return(res)
}

##' Modified Gene Set Enrichment Analysis (GSEA) of KEGG pathways.
##' The GSEA method is modified to accept geneList with large portion of zeros
##'
##' @title gseKEGG2
##' @param geneList order ranked geneList
##' @param organism supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
##' @param keyType one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
##' @param exponent weight of each step
##' @param nproc If not equal to zero, sets BPPARAM to use nproc workers (default = 1)
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @export
##' @return feaResult object
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



