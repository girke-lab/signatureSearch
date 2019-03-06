##' KEGG Enrichment Analysis of a gene set.
##' 
##' Given a vector of genes, this function will return the enrichment KEGG
##' categories with FDR control.
##'
##' @param gene a vector of entrez gene id.
##' @param organism supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
##' @param keyType one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param maxGSSize maximal size of genes annotated for testing
##' @param qvalueCutoff qvalue cutoff
##' @return A \code{feaResult} instance.
##' @importMethodsFrom AnnotationDbi mappedkeys
##' @importMethodsFrom AnnotationDbi mget
##' @importClassesFrom methods data.frame
enrichKEGG <- function(gene,
                       organism          = "hsa",
                       keyType           = "kegg",
                       pvalueCutoff      = 0.05,
                       pAdjustMethod     = "BH",
                       universe,
                       minGSSize         = 5,
                       maxGSSize         = 500,
                       qvalueCutoff      = 0.2) {
  
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
  res@ontology <- "KEGG"
  res@organism <- species
  return(res)
}

##' @importFrom clusterProfiler download_KEGG 
prepare_KEGG <- function(species, KEGG_Type="KEGG", keyType="kegg") {
  kegg <- download_KEGG(species, KEGG_Type, keyType)
  build_Anno(kegg$KEGGPATHID2EXTID,
             kegg$KEGGPATHID2NAME)
}

organismMapper <- function(organism) {
  ## it only map those previous supported organism
  
  if (organism == "anopheles") {
    species <- "aga"
  } else if (organism == "arabidopsis") {
    species <- "ath"
  } else if (organism == "bovine") {
    species <- "bta"
  } else if (organism == "canine") {
    species <- "cfa"
  } else if (organism == "chicken") {
    species <- "gga"
  } else if (organism == "chipm") {
    species <- "ptr"
  } else if (organism == "ecolik12") {
    species <- "eco"
  } else if (organism == "ecsakai") {
    species <- "ecs"
  } else if (organism == "fly") {
    species <- "dme"
  } else if (organism == "human") {
    species <- "hsa"
  } else if (organism == "malaria") {
    species <- "pfa"
  } else if (organism == "mouse") {
    species <- "mmu"
  } else if (organism == "pig") {
    species <- "ssc"
  } else if (organism == "rat") {
    species <- "rno"
  } else if (organism == "rhesus") {
    species <- "mcc"
  } else if (organism == "worm" || organism == "celegans") {
    species <- "cel"
  } else if (organism == "xenopus") {
    species <- "xla"
  } else if (organism == "yeast") {
    species <- "sce"
  } else if (organism == "zebrafish") {
    species <- "dre"
  } else {
    species <- organism
  }
  return(species)
}

