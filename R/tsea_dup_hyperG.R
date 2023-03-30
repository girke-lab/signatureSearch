#' @title FEA Methods
#' @rdname fea
#' @description 
#' The Target Set Enrichment Analysis (TSEA) with hypergeometric test 
#' (\code{tsea_dup_hyperG} function) performs TSEA based on a modified 
#' hypergeometric test that supports test sets with duplications. This is 
#' achieved by maintaining the frequency information of
#' duplicated items in form of weighting values. 
#' @details 
#' The classical hypergeometric test assumes uniqueness in its test sets. To
#' maintain the duplication information in the test sets used for TSEA, the 
#' values of the total number of genes/proteins in the test set and the number 
#' of genes/proteins in the test set annotated at a functional category are
#' adjusted by maintaining their frequency information in the test set rather 
#' than counting each entry only once. Removing duplications in TSEA would be 
#' inappropriate since it would erase one of the most important pieces of 
#' information of this approach.
#' @section Column description:
#' Descriptions of the columns in FEA result tables stored in the
#' \code{feaResult} object that can be accessed with the \code{result} method in
#' tabular format, here \code{tibble}.
#' \itemize{
#'     \item ont: in case of GO, one of BP, MF, CC, or ALL
#'     \item ID: GO or KEGG IDs
#'     \item Description: description of functional category
#'     \item GeneRatio: ratio of genes in the test set that are annotated at a 
#'     specific GO node or KEGG pathway
#'     \item BgRatio: ratio of background genes that are annotated
#'     at a specific GO node or KEGG pathway
#'     \item itemID: IDs of items (genes for TSEA, drugs for DSEA) overlapping 
#'     among test and annotation sets.
#'     \item setSize: size of the functional category
#'     \item pvalue from \code{tsea_dup_hyperG}: raw p-value of enrichment test
#'     \item p.adjust: p-value adjusted for multiple hypothesis testing based 
#'     on method specified under pAdjustMethod argument
#'     \item qvalue: q value calculated with R's qvalue function to control FDR
#'     \item enrichmentScore: ES from the GSEA algorithm 
#'     (Subramanian et al., 2005). The score is calculated by walking down the 
#'     gene list L, increasing a running-sum statistic when we encounter a gene 
#'     in S and decreasing when it is not. The magnitude of the increment 
#'     depends on the gene scores. The ES is the maximum deviation from zero 
#'     encountered in the random walk. It corresponds to a weighted 
#'     Kolmogorov-Smirnov-like statistic.
#'     \item NES: Normalized enrichment score. The positive and negative 
#'     enrichment scores are normalized separately by permutating the 
#'     composition of the gene list L \code{nPerm} times, and dividing the 
#'     enrichment score by the mean of the permutation ES with the same sign.
#'     \item pvalue from \code{tsea_mGSEA}: The nominal p-value of the ES is 
#'     calculated using a permutation test. Specifically, the composition of 
#'     the gene list L is permuted and the ES of the gene set is recomputed for 
#'     the permutated data generating a null distribution for the ES. 
#'     The p-value of the observed ES is then calculated relative to this 
#'     null distribution.
#'     \item leadingEdge: Genes in the gene set S (functional category) that 
#'     appear in the ranked list L at, or before, the point where the running 
#'     sum reaches its maximum deviation from zero. It can be interpreted as 
#'     the core of a gene set that accounts for the enrichment signal.
#'     \item ledge_rank: Ranks of genes in 'leadingEdge' in gene list L.
#'     \item mabs: given a scored ranked gene list \eqn{L}, \eqn{mabs(S)}
#'     represents the mean absolute scores of the genes in set \eqn{S}.
#'     \item Nmabs: normalized \eqn{mabs(S)}
#' }
#' 
#' @param drugs character vector containing drug identifiers used for functional
#' enrichment testing. This can be the top ranking drugs from a GESS result. 
#' Internally, drug test sets are translated to the corresponding target protein
#' test sets based on the drug-target annotations provided under the 
#' \code{dt_anno} argument.
#' @param universe character vector defining the universe of genes/proteins. If
#' set as 'Default', it uses all genes/proteins present in the corresponding
#' annotation system (e.g. GO, KEGG or Reactome). If 'type' is 'GO', it can be assigned
#' a custom vector of gene SYMBOL IDs. If 'type' is 'KEGG' or 'Reactome', the 
#' vector needs to contain Entrez gene IDs.
#' @param type one of `GO`, `KEGG` or `Reactome` if TSEA methods. \code{type}
#' can also be set as `MOA` is DSEA methods are used.
#' @param ont character(1). If type is `GO`, assign \code{ont} (ontology) one of
#' `BP`,`MF`, `CC` or `ALL`. If type is `KEGG` or `Reactome`, \code{ont} is ignored.
#' @param pAdjustMethod p-value adjustment method, 
#' one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'
#' @param pvalueCutoff double, p-value cutoff to return only enrichment results
#' for functional categories meeting a user definable confidence threshold
#' @param qvalueCutoff double, qvalue cutoff, similar to \code{pvalueCutoff}
#' @param minGSSize integer, minimum size of each gene set in annotation system.
#' Annotation categories with less than \code{minGSSize} genes/drugs will 
#' be ignored by enrichment test. If \code{type} is 'MOA', it may be beneficial
#' to set \code{minGSSize} to lower values (e.g. 2) than for other 
#' functional annotation systems. This is because certain MOA categories 
#' contain only 2 drugs.
#' @param maxGSSize integer, maximum size of each gene set in annotation system.
#' Annotation categories with more genes/drugs annotated than \code{maxGSSize} 
#' will be ignored by enrichment test.
#' @param dt_anno drug-target annotation source. It is the same argument as the
#' \code{database} argument of the \code{\link{get_targets}} function.
#' Usually, it is recommended to set the 'dt_anno' to 'all' since it provides 
#' the most complete drug-target annotations. Choosing a single
#' annotation source results in sparser drug-target annotations
#' (particularly CLUE), and thus less complete enrichment results.
#' @param readable TRUE or FALSE, it applies when type is `KEGG` or `Reactome`
#' indicating whether to convert gene Entrez ids to gene Symbols in the 'itemID' 
#' column in the result table.
#' @param nPerm integer defining the number of permutation iterations for 
#' calculating p-values
#' @param drugList named numeric vector, where the names represent drug labels 
#' and the numeric component scores. This can be all drugs of a GESS result that
#' are ranked by GESS scores, such as NCS scores from the LINCS method. 
#' Note, drugs with scores of zero are ignored by this method.
#' @param exponent integer value used as exponent in GSEA algorithm. It defines
#' the weight of the items in the item set \emph{S}. Note, in DSEA the items 
#' are drug labels, while it is gene labels in the original GSEA.
#' @return \code{\link{feaResult}} object, the result table contains the
#' enriched functional categories (e.g. GO terms or KEGG pathways) ranked by 
#' the corresponding enrichment statistic.
#' @seealso \code{\link{feaResult}}, 
#'          \code{\link[signatureSearchData]{GO_DATA_drug}}
#' @references 
#' GSEA algorithm: 
#' Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., 
#' Gillette, M. A., Mesirov, J. P. (2005). Gene set enrichment analysis: a 
#' knowledge-based approach for interpreting genome-wide expression profiles. 
#' Proceedings of the National Academy of Sciences of the United States of
#' America, 102(43), 15545-15550. URL: https://doi.org/10.1073/pnas.0506580102
#' 
#' MeanAbs algorithm: 
#' Fang, Z., Tian, W., & Ji, H. (2012). A network-based 
#' gene-weighting approach for pathway analysis. Cell Research, 22(3), 
#' 565-580. URL: https://doi.org/10.1038/cr.2011.149
#' @import org.Hs.eg.db
#' @examples 
#' 
#' ############### TSEA dup_hyperG method ########
#' ## GO annotation system
#' # res1 <- tsea_dup_hyperG(drugs=drugs10, universe="Default", 
#' #                         type="GO", ont="MF", pvalueCutoff=0.05,
#' #                         pAdjustMethod="BH", qvalueCutoff=0.1, 
#' #                         minGSSize=5, maxGSSize=500)
#' # result(res1)
#' #
#' ## KEGG annotation system
#' # res2 <- tsea_dup_hyperG(drugs=drugs10, type="KEGG", 
#' #                         pvalueCutoff=0.1, qvalueCutoff=0.2, 
#' #                         minGSSize=10, maxGSSize=500)
#' #
#' ## Reactome annotation system
#' # res3 <- tsea_dup_hyperG(drugs=drugs10, type="Reactome", 
#' #                         pvalueCutoff=1, qvalueCutoff=1)
#' @export
tsea_dup_hyperG <- function(drugs, universe="Default", 
                            type="GO", ont="MF", 
                            pAdjustMethod="BH", pvalueCutoff=0.05, 
                            qvalueCutoff=0.05, 
                            minGSSize=5, maxGSSize=500,
                            dt_anno="all", readable=FALSE){
  # message("The query drugs [", length(drugs),"] are: \n", 
  #         paste0(drugs[seq_len(min(length(drugs), 10))], sep ="  "), "...")

  if(!any(type %in% c("GO", "KEGG", "Reactome"))){
    stop('"type" argument needs to be one of "GO", "KEGG" or "Reactome"')
  }
  drugs <- unique(tolower(drugs))
  targets <- get_targets(drugs, database = dt_anno)
  gnset <- na.omit(unlist(lapply(targets$t_gn_sym, function(i) 
    unlist(strsplit(as.character(i), split = "; ")))))
  
  OrgDb <- load_OrgDb("org.Hs.eg.db")
  valid_sym <- keys(OrgDb, keytype="SYMBOL")
  if(sum(gnset %in% valid_sym) == 0){
    return(NULL)
  }
  
  if(type=="GO"){
    if(universe=="Default"){
      # Get universe genes in GO annotation system
      universe <- univ_go
    }
    ego <- enrichGO2(gene = gnset, universe = universe, OrgDb = "org.Hs.eg.db", 
                    keytype = 'SYMBOL', ont = ont, 
                    pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, 
                    qvalueCutoff = qvalueCutoff, 
                    minGSSize = minGSSize, maxGSSize = maxGSSize)
    if(! is.null(ego)){
        drugs(ego) <- drugs
    }
    return(ego)
  }
  
  gnset_entrez <- na.omit(suppressMessages(
    AnnotationDbi::select(OrgDb, keys = gnset, keytype = "SYMBOL", 
                          columns = "ENTREZID")$ENTREZID))
  
  if(type=="KEGG"){
    if(universe=="Default"){
      KEGG_DATA <- prepare_KEGG(species="hsa", "KEGG", keyType="kegg")
      keggterms <- get("PATHID2EXTID", KEGG_DATA)
      universe <- unique(unlist(keggterms))
    }
    eres <- enrichKEGG2(gene=gnset_entrez, organism="hsa", keyType="kegg", 
                      pvalueCutoff=pvalueCutoff, qvalueCutoff=qvalueCutoff,
                      pAdjustMethod=pAdjustMethod, universe=universe, 
                      minGSSize=minGSSize, maxGSSize=maxGSSize, readable=readable)
  }
  
  if(type=="Reactome"){
    if(universe=="Default"){
      Reactome_DATA <- get_Reactome_DATA(organism="human")
      raterms <- get("PATHID2EXTID", Reactome_DATA)
      universe <- unique(unlist(raterms))
    }

    eres <- enrichReactome(gene=gnset_entrez, organism="human",
                           pvalueCutoff=pvalueCutoff, qvalueCutoff=qvalueCutoff, 
                           pAdjustMethod=pAdjustMethod, universe=universe,
                           minGSSize=minGSSize, maxGSSize=maxGSSize, readable=readable)
  }

  if(!is.null(eres)){
    drugs(eres) <- drugs
  }
  return(eres)
}
