##' This modified Gene Set Enrichment Analysis (GSEA) of GO terms supports
##' gene test sets with large numbers of zeros.
##'
##' @title Modified GSEA with GO Terms
##' @param geneList named numeric vector with gene SYMBOLs in the name slot
##' decreasingly ranked by scores in the data slot.
##' @param ont one of "BP", "MF", "CC" or "ALL"
##' @param OrgDb OrgDb, e.g., "org.Hs.eg.db".
##' @param keyType keytype of gene
##' @param exponent weight of each step
##' @param nproc if not equal to zero, sets \code{BPPARAM} to use \code{nproc} 
##' workers (default = 1)
##' @param nPerm permutation numbers
##' @param minGSSize integer, minimum size of each gene set in annotation system
##' @param maxGSSize integer, maximum size of each gene set in annotation system
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @return feaResult object
##' @examples 
##' data(targetList)
##' # gsego <- gseGO2(geneList=targetList, ont="MF", OrgDb="org.Hs.eg.db",
##' #                 pvalueCutoff=1)
##' # head(gsego)
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
    
    res <- GSEA_internal2(geneList = geneList,
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
    res <- select_ont(res, ont, GO_DATA)
    ont(res) <- ont
    og(res) <- get_organism(OrgDb)
    return(res)
}

##' This modified Gene Set Enrichment Analysis (GSEA) of KEGG pathways supports
##' gene test sets with large numbers of zeros.
##'
##' @title Modified GSEA with KEGG
##' @param geneList named numeric vector with gene ids in the name slot 
##' decreasingly ranked by scores in the data slot.
##' @param organism supported organism listed in
##' URL: http://www.genome.jp/kegg/catalog/org_list.html
##' @param keyType one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
##' @param exponent weight of each step
##' @param nproc if not equal to zero, sets \code{BPPARAM} to use \code{nproc} 
##' workers (default = 1)
##' @param nPerm permutation numbers
##' @param minGSSize integer, minimum size of each gene set in annotation system
##' @param maxGSSize integer, maximum size of each gene set in annotation system
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @param readable TRUE or FALSE indicating whether to convert gene Entrez ids
##' to gene Symbols in the 'itemID' column in the FEA result table.
##' @return feaResult object
##' @examples 
##' # Gene Entrez id should be used for KEGG enrichment
##' data(geneList, package="DOSE")
##' #geneList[100:length(geneList)]=0
##' #gsekk <- gseKEGG2(geneList=geneList, pvalueCutoff = 1)
##' #head(gsekk)
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
                    verbose           = TRUE, readable=FALSE) {
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
    if(readable) result(res) <- set_readable(result(res), geneCol="leadingEdge")
    og(res) <- species
    ont(res) <- "KEGG"

    return(res)
}

##' This modified Gene Set Enrichment Analysis (GSEA) of Reactome pathways
##' supports gene test sets with large numbers of zeros.
##'
##' @title Modified GSEA with Reactome
##' @param geneList order ranked geneList
##' @param organism one of "human", "rat", "mouse", "celegans", "yeast", 
##' "zebrafish", "fly".
##' @param exponent integer value used as exponent in GSEA algorithm.
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of each geneSet for analyzing
##' @param nPerm integer defining the number of permutation iterations for 
##' calculating p-values
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' TRUE or FALSE indicating whether to convert gene Entrez ids
##' to gene Symbols in the 'itemID' column in the FEA result table.
##' @param readable TRUE or FALSE indicating whether to convert gene Entrez ids
##' to gene Symbols in the 'itemID' column in the FEA result table.
##' @return feaResult object
##' @examples 
##' # Gene Entrez id should be used for Reactome enrichment
##' data(geneList, package="DOSE")
##' #geneList[100:length(geneList)]=0
##' #rc <- gseReactome(geneList=geneList, pvalueCutoff=1)
##' @export
gseReactome <- function(geneList,
                        organism      = "human",
                        exponent      = 1,
                        nPerm         = 1000,
                        minGSSize     = 10,
                        maxGSSize     = 500,
                        pvalueCutoff  = 0.05,
                        pAdjustMethod = "BH",
                        verbose       = TRUE, readable=FALSE){
    Reactome_DATA <- get_Reactome_DATA(organism)
    
    res <- GSEA_internal2(geneList = geneList,
                          exponent = exponent,
                          nPerm = nPerm,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          verbose = verbose,
                          USER_DATA = Reactome_DATA,
                          nproc = 1)
    
    if (is.null(res))
        return(res)
    if(readable) result(res) <- set_readable(result(res), geneCol="leadingEdge")
    og(res) <- organism
    ont(res) <- "Reactome"
    return(res)
}
