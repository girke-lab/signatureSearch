##' A drug-target interaction network can be built if a drug set and a 
##' target/protein set is provided. If the target set is a GO category or a 
##' KEGG pathway, their term ID could be directly used.
##' @title plot drug-target interaction network
##' @param drugs A character vector representing drugs used to build a 
##' drug-target network 
##' @param set could be GO term ID, KEGG pathway ID or a character vector of 
##' gene set with SYMBOL ids.
##' @param ont if the `set` is a GO term ID, `ont` is the ontology that GO term 
##' is belong to. One of "BP", "MF", "CC" or "ALL"
##' @param ... additional parameters
##' @return visNetwork plot
##' @import visNetwork
##' @importFrom AnnotationDbi select
##' @importFrom scales cscale
##' @importFrom scales seq_gradient_pal
##' @examples 
##' data(drugs)
##' dtnetplot(drugs=drugs, set="GO:0032041", ont = "MF", 
##'     desc="NAD-dependent histone deacetylase activity (H3-K14 specific)")
##' @export dtnetplot
##' 
dtnetplot <- function(drugs, set, ont = NULL, ...) {
  if(grepl("GO:\\d{7}",set)){
    ont %<>% toupper
    if(is.null(ont) | !any(ont %in% c("BP","MF","CC","ALL"))) 
      stop("The 'set' is a GO term ID, please set 'ont' as one of 
           BP, MF, CC or ALL")
    # download goAnno.rds and save it to cache
    fl <- download_data_file(url=
    "http://biocluster.ucr.edu/~yduan004/signatureSearch_data/goAnno.rds",
                             rname="goAnno")
    goAnno <- readRDS(fl)
    
    go_gene <- unique(goAnno$SYMBOL[goAnno$ONTOLOGYALL == ont & 
                                      goAnno$GOALL == set])
  }
  if(grepl("hsa\\d{5}",set)){
    KEGG_DATA <- prepare_KEGG(species="hsa", "KEGG", keyType="kegg")
    p2e <- get("PATHID2EXTID", envir=KEGG_DATA)
    go_gene_entrez = p2e[[set]]
    # convert Entrez ids in KEGG pathways to gene SYMBOL
    db <- load_OrgDb("org.Hs.eg.db")
    go_gene_map <- suppressMessages(
      select(db, keys = go_gene_entrez, keytype = "ENTREZID", columns="SYMBOL"))
    go_gene <- unique(go_gene_map$SYMBOL)
  }
  
  # get drug targets in DrugBank, STITCH, LINCS and calculate targets weight
  dtslash <- get_targets(drugs, database="all")
  dtlink <- slash2link(dtslash)
  dtlink_go <- dtlink[dtlink$t_gn_sym %in% go_gene,]
  go_gene_nottar <- setdiff(go_gene, unique(dtlink$t_gn_sym))
  go_gene_tar <- intersect(go_gene, unique(dtlink$t_gn_sym))
  
  message(length(go_gene_tar), "/", length(go_gene), " (", 
          round(length(go_gene_tar)/length(go_gene)*100,2),
          "%) genes in the gene set are targeted by query drugs, which are ", 
          paste0(go_gene_tar, collapse = " / "), "\n")
  if(length(go_gene_nottar>0)){
    message(length(go_gene_nottar), "/", length(go_gene), " (", 
        round(length(go_gene_nottar)/length(go_gene)*100,2),
        "%) genes in the gene set are not targeted by query drugs, which are ", 
        paste0(go_gene_nottar, collapse = " / "), "\n")
  }
  
  drugs_tar <- unique(dtlink_go$drug_name)
  drugs_no <- setdiff(drugs, drugs_tar)
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
  visNetwork(nodes, edges, height = "500px", width = "100%") %>%
    visLegend(width = 0.05, position = "right", addNodes = lnodes, 
              useGroups = FALSE) %>%
    visOptions(highlightNearest=list(enabled=TRUE, degree=1, hover=TRUE), 
               nodesIdSelection = TRUE)
}
