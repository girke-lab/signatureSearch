##' Functional modules of GESS and FEA results can be rendered as interactive
##' drug-target networks using the \code{dtnetplot} function form
##' \code{signatureSearch}. For this, a character vector of drug names along 
##' with an identifier of a chosen functional category are passed on to the 
##' drugs and set arguments, respectively. The resulting plot depicts the 
##' corresponding drug-target interaction network. Its interactive features 
##' allow the user to zoom in and out of the network, and to select network 
##' components in the drop-down menu located in the upper left corner of the 
##' plot.
##' @title Drug-Target Network Visualization
##' @param drugs A character vector of drug names
##' @param set character(1) GO term ID or KEGG pathway ID. Alternatively, a
##' character vector of gene SYMBOLs can be assigned.
##' @param ont if `set` is a GO term ID, `ont` is the corresponding ontology 
##' that GO term belongs to. One of 'BP', 'MF' or 'CC'
##' @param desc character(1), description of the chosen functional category or 
##' target set
##' @param verbose TRUE or FALSE, whether to print messages
##' @param ... additional parameters of \code{\link[visNetwork]{visNetwork}} 
##' function.
##' @return visNetwork plot
##' @import visNetwork
##' @importFrom AnnotationDbi select
##' @importFrom scales cscale
##' @importFrom scales seq_gradient_pal
##' @examples 
##' data(drugs10)
##' dtnetplot(drugs=drugs10, 
##'     set=c("HDAC1", "HDAC2", "HDAC3", "HDAC11", "FOX2"),
##'     desc="NAD-dependent histone deacetylase activity (H3-K14 specific)")
##' @export dtnetplot
##' 
dtnetplot <- function(drugs, set, ont=NULL, desc=NULL, verbose=FALSE, ...) {
  if(grepl("GO:\\d{7}",set)[1]){
    ont %<>% toupper
    if(is.null(ont) | !any(ont %in% c("BP","MF","CC"))) 
      stop("The 'set' is a GO term ID, please set 'ont' as one of 
           BP, MF or CC")
    # download goAnno.rds and save it to cache
    goAnno <- suppressMessages(ah[["AH69084"]])
    go_gene <- unique(goAnno$SYMBOL[goAnno$ONTOLOGYALL == ont & 
                                      goAnno$GOALL == set])
  } else if(grepl("hsa\\d{5}",set)[1]){
    KEGG_DATA <- prepare_KEGG(species="hsa", "KEGG", keyType="kegg")
    p2e <- get("PATHID2EXTID", envir=KEGG_DATA)
    go_gene_entrez = p2e[[set]]
    # convert Entrez ids in KEGG pathways to gene SYMBOL
    go_gene_map <- suppressMessages(
      select(org.Hs.eg.db, keys = go_gene_entrez, keytype = "ENTREZID", 
             columns="SYMBOL"))
    go_gene <- unique(go_gene_map$SYMBOL)
  } else {
      go_gene <- set
  }
  
  # get drug targets in DrugBank, STITCH, LINCS and calculate targets weight
  dtslash <- get_targets(drugs, database="all")
  dtlink <- slash2link(dtslash)
  dtlink_go <- dtlink[dtlink$t_gn_sym %in% go_gene,]
  go_gene_nottar <- setdiff(go_gene, unique(dtlink$t_gn_sym))
  go_gene_tar <- intersect(go_gene, unique(dtlink$t_gn_sym))
  if(verbose){
    message(length(go_gene_tar), "/", length(go_gene), " (", 
            round(length(go_gene_tar)/length(go_gene)*100,2),
            "%) genes in the gene set are targeted by query drugs, which are ", 
            paste0(go_gene_tar, collapse = " / "), "\n")
    if(length(go_gene_nottar>0)){
      message(length(go_gene_nottar), "/", length(go_gene), " (", 
              round(length(go_gene_nottar)/length(go_gene)*100,2),
              "%) genes in the gene set are not targeted by query drugs, ",
              "which are ", paste0(go_gene_nottar, collapse = " / "), "\n")
      }
  }
  
  drugs_tar <- unique(dtlink_go$drug_name)
  drugs_no <- setdiff(drugs, drugs_tar)
  if(verbose){
    message(length(drugs_tar), "/", length(drugs), " (", 
            round(length(drugs_tar)/length(drugs)*100,2),
            "%) drugs target genes/proteins in the gene set. They are ", 
            paste0(drugs_tar, collapse = " / "), "\n")
    if(length(drugs_no>0)){
      message(length(drugs_no), "/", length(drugs), " (", 
              round(length(drugs_no)/length(drugs)*100,2),
              "%) drugs don't target genes/proteins in the gene set. They are ",
              paste0(drugs_no, collapse = " / "), "\n")
    }
  }
  
  ## scale node colors based on targetWeight and drugWeight(dw)
  tw <- table(dtlink_go$t_gn_sym)
  tw_vec <- as.numeric(tw); names(tw_vec) <- names(tw)
  tw_vec0 <- rep(0, length(go_gene_nottar))
  names(tw_vec0) <- go_gene_nottar
  targetWeight <- c(tw_vec, tw_vec0)
  col.scale_tar <- cscale(targetWeight, seq_gradient_pal("#B3B3B3","#FF5C32"))
  
  dw <- table(dtlink_go$drug_name)
  dw_vec <- as.numeric(dw); names(dw_vec) <- names(dw)
  col.scale_drug <- cscale(dw_vec, seq_gradient_pal("#ffd07aff", "orange"))
  
  # Use visNetwork to plot
  lnodes <- data.frame(label = c("drugs", "targets"),
                       shape = c("box", "circle"), 
                       color = c("orange", "#FF5C32"),
                       title = "Groups", id = seq_len(2))
  nodes <- data.frame(id = c(drugs_tar, go_gene), label = c(drugs_tar, go_gene),
                      group = c(rep("drugs", length(drugs_tar)), 
                                rep("targets", length(go_gene))),
                      shape = c(rep("box", length(drugs_tar)), 
                                rep("circle", length(go_gene))), 
                      color = c(col.scale_drug[drugs_tar], 
                                col.scale_tar[go_gene]),
                      value = c(5*table(dtlink_go$drug_name), 
                                rep(5, length(go_gene))))
  edges <- dtlink_go; colnames(edges) = c("from", "to")
  visNetwork(nodes, edges, width = "100%", ...) %>%
    visLegend(width = 0.05, position = "right", addNodes = lnodes, 
              useGroups = FALSE) %>%
    visOptions(highlightNearest=list(enabled=TRUE, degree=1, hover=TRUE), 
               nodesIdSelection = TRUE)
}
