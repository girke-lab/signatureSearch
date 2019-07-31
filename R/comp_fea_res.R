##' Plot top re-ranking functional categories from different enrichment methods.
##' 
##' The `comp_fea_res` function re-ranks the functional categories from 
##' differentFEA methods by using the mean rank of each functional category
##' across the 5 FEA methods, the lower of the mean rank, the topper of the 
##' functional category. Since the functional categories have different number 
##' of support methods, the mean rank of a functional category is adjusted by 
##' multiplying the total number of methods and divided by the number of support 
##' methods. In this way, if a functional category is supported by only one 
##' method, the mean rank will be increased 5 times for penalty. 
##'
##' @title Plot Comparing FEA Results
##' @param table_list a named list of tables extracted from feaResults from 
##' different FEA methods. The names of the list is required, which could be 
##' name of the enrichment methods.  
##' @param rank_stat character(1), column name of the enrichment statisic used 
##' for ranking the functional categories, e.g. 'pvalue' or 'p.adjust'. Note,
##' the column name should exist in each table in the 'table_list'.
##' @param Nshow integer, number of top re-ranking functional categories shown 
##' in the plot
##' @importFrom ggplot2 aes_
##' @importFrom DOSE theme_dose
##' @return igraph object
##' @examples 
##' data(drugs)
##' # table_list = list("dup_hyperG" = result(dup_hyperG_res), 
##' #                   "mGSEA" = result(mgsea_res), 
##' #                   "mabs" = result(mabs_res), 
##' #                   "hyperG" = result(hyperG_res), 
##' #                   "GSEA" = result(gsea_res))
##' # comp_fea_res(table_list, rank_stat="pvalue", Nshow=20)
##' @export

comp_fea_res <- function(table_list, rank_stat="pvalue", Nshow=20){
    if(is.null(names(table_list))){
        stop(paste('The "table_list" should be a list with names,', 
             'which could be names of the enrichment methods'))
    }
    df <- NULL
    for(i in seq_along(table_list)){
        tb <- table_list[[i]]
        tb_order <- tb[order(tb[[rank_stat]]),]
        tb2 <- data.frame(tb_order[, c("Description", rank_stat)], 
                             method=names(table_list)[i], 
                             rank=seq_len(nrow(tb_order)))
        df <- rbind(df, tb2)
    }
    # re-ranking according to mean ranks and adjusted for number of support
    # methods (*5/number of support methods)
    cat_ranks_list <- split(df$rank, df$Description) 
    cat_mrk_adj <- sapply(cat_ranks_list, function(i){
        m <- mean(i)*5/length(i)
    })
    cat_rerk <- cat_mrk_adj[order(cat_mrk_adj)]
    cat_top <- names(cat_rerk)[seq_len(Nshow)]
    df2 <- df[df$Description %in% cat_top,]
    df2$Description <- factor(df2$Description, levels = rev(cat_top), 
                            ordered = TRUE)
    p <- ggplot(df2, aes_(x = ~method, y = ~Description)) +
    geom_point(size=4) + aes_string(color=rank_stat) + 
    scale_colour_gradient(low="red", high="blue") +
    theme_dose(font.size=12) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
    return(p)
}


