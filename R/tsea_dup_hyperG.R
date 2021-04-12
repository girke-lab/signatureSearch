#' Target Set Enrichment Analysis (TSEA) with Hypergeometric Test
#' 
#' The \code{tsea_dup_hyperG} function performs Target Set Enrichment Analysis
#' (TSEA) based on a modified hypergeometric test that supports test sets with
#' duplications. This is achieved by maintaining the frequency information of
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
#' The TSEA results (including \code{tsea_dup_hyperG}) stored in the
#' \code{feaResult} object can be returned with the \code{result} method in
#' tabular format, here \code{tibble}. The columns of this \code{tibble} are
#' described below.
#' \itemize{
#'     \item GeneRatio: ratio of genes in the test set that are annotated at a 
#'     specific GO node or KEGG pathway
#'     \item BgRatio: ratio of background genes that are annotated
#'     at a specific GO node or KEGG pathway
#'     \item pvalue: raw p-value of enrichment test
#' }
#' Additional columns are described under the 'result' slot of the
#' \code{\link{feaResult}} object.
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
#' @param type one of `GO`, `KEGG` or `Reactome`
#' @param ont character(1). If type is `GO`, assign \code{ont} (ontology) one of
#' `BP`,`MF`, `CC` or `ALL`. If type is `KEGG` or `Reactome`, \code{ont} is ignored.
#' @param pAdjustMethod p-value adjustment method, 
#' one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'
#' @param pvalueCutoff double, p-value cutoff
#' @param qvalueCutoff double, qvalue cutoff
#' @param minGSSize integer, minimum size of each gene set in annotation system
#' @param maxGSSize integer, maximum size of each gene set in annotation system
#' @param dt_anno drug-target annotation source. Currently, one of 'DrugBank',
#' 'CLUE', 'STITCH' or 'all'. If 'dt_anno' is 'all', the targets from the
#' DrugBank, CLUE and STITCH databases will be combined. Usually, it is
#' recommended to set the 'dt_anno' to 'all' since it provides the most
#' complete drug-target annotations. Choosing a single
#' annotation source results in sparser drug-target annotations
#' (particularly CLUE), and thus less complete enrichment results.
#' @param readable TRUE or FALSE, it applies when type is `KEGG` or `Reactome`
#' indicating whether to convert gene Entrez ids to gene Symbols in the 'itemID' 
#' column in the result table.
#' @return \code{\link{feaResult}} object, the result table contains the
#' enriched functional categories (e.g. GO terms or KEGG pathways) ranked by 
#' the corresponding enrichment statistic.
#' @seealso \code{\link{feaResult}}, \code{\link{fea}}
#' @examples 
#' data(drugs10)
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
