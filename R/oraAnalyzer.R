##' Given a vector of gene identifiers, this function returns GO term enrichment
##' results based on a hypergeometric test with duplication support in the test 
##' set.
##' 
##' @title GO Term Enrichment with Hypergeometric Test
##' @param gene a vector of gene SYMBOL ids (here the test set)
##' @param OrgDb OrgDb
##' @param keytype Gene ID type of test set
##' @param ont One of "MF", "BP", "CC" or "ALL"
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", 
##' "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param qvalueCutoff qvalue cutoff
##' @param minGSSize minimum size of each gene set in annotation system
##' @param maxGSSize maximum size of each gene set in annotation system
##' @param pool If ont='ALL', whether 3 GO ontology should be combined
##' @return A \code{feaResult} instance.
##' @seealso \code{\link{feaResult-class}}
##' @examples 
##' # The method supports duplicated elements in 'gene', 
##' # which should be gene SYMBOL ids for GO term enrichment.
##' gene <- c(rep("HDAC1",4), rep("HDAC3",2), "SOX8", "KLK14")
##' # data(targetList)
##' # ego <- enrichGO2(gene = gene, OrgDb="org.Hs.eg.db", ont="MF",
##' #                  universe=names(targetList))
##' @export
enrichGO2 <- function(gene,
                      OrgDb,
                      keytype = "SYMBOL",
                      ont="MF",
                      pvalueCutoff=0.05,
                      pAdjustMethod="BH",
                      universe,
                      qvalueCutoff = 0.2,
                      minGSSize = 5,
                      maxGSSize = 500,
                      pool=FALSE) {
    
    ont %<>% toupper
    ont <- match.arg(ont, c("BP", "CC", "MF", "ALL"))
    # GO_DATA <- clusterProfiler:::get_GO_data(OrgDb, ont, keytype)
    # download GO_DATA.rds from AnnotationHub to save time by avoiding 
    # building GO_DATA from scratch
    eh <- suppressMessages(ExperimentHub())
    GO_DATA <- suppressMessages(eh[["EH3231"]])
    if (missing(universe))
        universe <- NULL
    
    res <- enricher_internal(gene,
                             pvalueCutoff=pvalueCutoff,
                             pAdjustMethod=pAdjustMethod,
                             universe = universe,
                             qvalueCutoff = qvalueCutoff,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             USER_DATA = GO_DATA)
    
    if (is.null(res))
        return(res)
    res <- select_ont(res, ont, GO_DATA)
    og(res) <- get_organism(OrgDb)
    ont(res) <- ont
    return(res)
}

##' Given a vector of gene identifiers, this function returns KEGG pathway 
##' enrichment results based on a hypergeometric test with duplication support 
##' in the test set.
##' 
##' @title KEGG Pathway Enrichment with Hypergeometric Test
##' @param gene a vector of entrez gene ids (here the test set)
##' @param organism supported organism are listed in 
##' http://www.genome.jp/kegg/catalog/org_list.html
##' @param keyType one of "kegg", 'ncbi-geneid', 'ncbi-proteinid' or 'uniprot'
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", 
##' "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by 
##' ontology term for testing.
##' @param maxGSSize maximal size of genes annotated for testing
##' @param qvalueCutoff qvalue cutoff
##' @param readable TRUE or FALSE indicating whether to convert gene Entrez ids
##' to gene Symbols in the 'itemID' column in the FEA result table.
##' @return A \code{feaResult} instance.
##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importMethodsFrom AnnotationDbi mget
##' @importClassesFrom methods data.frame
##' @examples 
##' # Method supports duplicated elements in "gene", which should be entrez ids
##' gene <- c(rep("4312",4), rep("8318",2), "991", "10874")
##' #data(geneList, package="DOSE")
##' #kk <- enrichKEGG2(gene = gene, universe=names(geneList))
##' #head(kk)
##' @export
enrichKEGG2 <- function(gene,
                        organism          = "hsa",
                        keyType           = "kegg",
                        pvalueCutoff      = 0.05,
                        pAdjustMethod     = "BH",
                        universe,
                        minGSSize         = 5,
                        maxGSSize         = 500,
                        qvalueCutoff      = 0.2, readable=FALSE) {
    
    species <- organismMapper(organism)
    KEGG_DATA <- prepare_KEGG(species, "KEGG", keyType)
    res <- enricher_internal(gene,
                             pvalueCutoff  = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             universe      = universe,
                             minGSSize     = minGSSize,
                             maxGSSize     = maxGSSize,
                             qvalueCutoff  = qvalueCutoff,
                             USER_DATA = KEGG_DATA)
    if (is.null(res))
        return(res)
    if(readable) result(res) <- set_readable(result(res))
    ont(res) <- "KEGG"
    og(res) <- species
    return(res)
}

##' Given a vector of gene identifiers, this function returns MOA category 
##' enrichment results based on a hypergeometric test with duplication support 
##' in the test set. The universe for the test is set to the unique genes 
##' encoding the target proteins present in the MOA annotation system from the 
##' ChEMBL database. 
##' 
##' @title MOA Category Enrichment with Hypergeometric Test
##' @param gene a vector of entrez gene ids (here the test set)
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", 
##' "bonferroni", "BH", "BY", "fdr", "none"
##' @param qvalueCutoff qvalue cutoff
##' @return A \code{feaResult} instance.
##' @seealso \code{\link{feaResult-class}}
##' @examples 
##' data(geneList, package="DOSE")
##' emoa <- enrichMOA(gene = names(geneList)[seq(3)])
##' head(emoa)
##' @export
enrichMOA <- function(gene,
                      pvalueCutoff=0.05,
                      pAdjustMethod="BH",
                      qvalueCutoff = 0.2) {
    data("chembl_moa_list", envir = environment())
    MOA_DATA_chembl <- get_MOA_data(chembl_moa_list, keytype="entrez")
    # get all the gene entrez ids in the MOA annotation system as universe
    ext2path <- get("EXTID2PATHID", envir = MOA_DATA_chembl)
    universe <- names(ext2path)
    
    res <- enricher_internal(gene,
                             pvalueCutoff=pvalueCutoff,
                             pAdjustMethod=pAdjustMethod,
                             universe = universe,
                             qvalueCutoff = qvalueCutoff,
                             minGSSize = 0,
                             maxGSSize = 100,
                             USER_DATA = MOA_DATA_chembl)
    
    if (is.null(res))
        return(res)
    og(res) <- "Homo sapiens"
    ont(res) <- "MOA"
    return(res)
}

##' Reactome Enrichment Analysis of a gene set.
##' Given a vector of genes, this function will return the enriched Reactome 
##' pathways with FDR control from hypergeometric test.
##'
##' @param gene a vector of entrez gene id.
##' @param organism one of "human", "rat", "mouse", "celegans", "yeast", 
##' "zebrafish", "fly".
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", 
##' "BH", "BY", "fdr", "none"
##' @param qvalueCutoff Cutoff value of qvalue
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by functional term for testing.
##' @param maxGSSize maximal size of each gene set for analyzing
##' @param readable TRUE or FALSE indicating whether to convert gene Entrez ids
##' to gene Symbols in the 'itemID' column in the FEA result table.
##' @return A \code{feaResult} instance.
##' @export
##' @seealso \code{\link{feaResult-class}}
##' @examples
##' # This method supports duplicated elements in "gene"
##' gene <- c(rep("4312",4), rep("8318",2), "991", "10874")
##' #data(geneList, package="DOSE")
##' #rc <- enrichReactome(gene=gene, universe=names(geneList))
##' #result(rc)
enrichReactome <- function(gene, organism="human",
                           pvalueCutoff=0.05, pAdjustMethod="BH",
                           qvalueCutoff=0.2, universe,
                           minGSSize=5, maxGSSize=500, readable=FALSE){
    
    Reactome_DATA <- get_Reactome_DATA(organism)
    
    res <- enricher_internal(gene,
                             pvalueCutoff=pvalueCutoff,
                             pAdjustMethod=pAdjustMethod,
                             qvalueCutoff=qvalueCutoff,
                             universe = universe,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             USER_DATA = Reactome_DATA)
    
    if (is.null(res))
        return(res)
    if(readable) result(res) <- set_readable(result(res))
    og(res) <- organism
    ont(res) <- "Reactome"
    return(res)
}
