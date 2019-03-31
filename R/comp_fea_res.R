##' Plot consistency comparison of enriched terms from FEA methods
##' 
##' The enriched GO/KEGG categories were re-ranked by putting the terms 
##' supported by more methods to topper. For terms supported by the same number
##' of methods, the lower of the mean ranks of the terms across 
##' enrichment methods, the topper of the terms. 
##'
##' @title plot comparison of FEA enrichment results
##' @param dup_hyperG_res feaResult object from 'dup_hyperG' method
##' @param mgsea_res feaResult object from 'mGSEA' method
##' @param mabs_res feaResult object from 'mabs' method
##' @param hyperG_res feaResult object from 'hyperG' method
##' @param gsea_res feaResult object from 'GSEA' method
##' @param Ntop number of top enriched categories from each method used to do 
##' comparison. By default, the plot only compare the top 20 ranking terms
##' @param type "GO" or "KEGG", indicating the type of enriched terms
##' @importFrom ggplot2 aes_
##' @importFrom DOSE theme_dose
##' @return igraph object
##' @examples 
##' data(drugs)
##' #comp_fea_res(dup_hyperG_res, mgsea_res, mabs_res, 
##' #             hyperG_res, gsea_res, Ntop=20, type="GO")
##' @export

comp_fea_res <- function(dup_hyperG_res=NULL, mgsea_res=NULL, mabs_res=NULL, 
                         hyperG_res=NULL, gsea_res=NULL, Ntop=20, type="GO"){
  df <- NULL
  if(! is.null(dup_hyperG_res)){
    dup_tb <- result(dup_hyperG_res)
    dup_df <- data.frame(dup_tb[seq_len(min(nrow(dup_tb), Ntop)),
                                c("Description","p.adjust")], 
                         method="dup_hyperG", 
                         rank=seq_len(min(nrow(dup_tb), Ntop)))
    df <- rbind(df, dup_df)
  }
  if(! is.null(mgsea_res)){
    mgsea_tb <- result(mgsea_res)
    mgsea_df <- data.frame(mgsea_tb[seq_len(min(nrow(mgsea_tb), Ntop)),
                                    c("Description","p.adjust")], 
                           method="mGSEA", 
                           rank=seq_len(min(nrow(mgsea_tb), Ntop)))
    df <- rbind(df, mgsea_df)
  }
  if(! is.null(mabs_res)){
    mabs_tb <- result(mabs_res)
    mabs_df <- data.frame(mabs_tb[seq_len(min(nrow(mabs_tb), Ntop)),
                                  c("Description","p.adjust")], 
                          method="mabs", 
                          rank=seq_len(min(nrow(mabs_tb), Ntop)))
    df <- rbind(df, mabs_df)
  }
  if(! is.null(hyperG_res)){
    hpg_tb <- result(hyperG_res)
    hpg_df <- data.frame(hpg_tb[seq_len(min(nrow(hpg_tb), Ntop)),
                                c("Description","p.adjust")], 
                         method="hyperG", 
                         rank=seq_len(min(nrow(hpg_tb), Ntop)))
    df <- rbind(df, hpg_df)
  }
  if(! is.null(gsea_res)){
    gsea_tb <- result(gsea_res)
    gsea_df <- data.frame(gsea_tb[seq_len(min(nrow(gsea_tb), Ntop)),
                                  c("Description","p.adjust")], 
                          method="GSEA", 
                          rank=seq_len(min(nrow(gsea_tb), Ntop)))
    df <- rbind(df, gsea_df)
  }
  # rank descriptions
  ## first rank according to frequence
  freq <- as.factor(sort(table(df$Description), 
                         decreasing = TRUE)[seq_len(Ntop)])
  ## In each frequence, rank by mean of ranks
  desc_order <- NULL
  for (i in rev(levels(freq))){
    sub_df <- df[df$Description %in% names(freq[freq==i]),]
    sub_desc <- sort(tapply(sub_df$rank, sub_df$Description, mean))
    desc_order <- c(desc_order, names(sub_desc))
  }
  df2 <- df[df$Description %in% desc_order,]
  df2$Description <- factor(df2$Description, levels = rev(desc_order), 
                            ordered = TRUE)
  title=paste0("Consistency comparison of top ", Ntop, 
               " enriched \n ", type," terms from 5 methods")
  p <- ggplot(df2, aes_(x = ~method, y = ~Description)) +
    geom_point(size=4) + aes_string(color="p.adjust") + 
    scale_colour_gradient(low="red", high="blue") +
    ggtitle(title) + theme_dose(font.size=12)
  return(p)
}


