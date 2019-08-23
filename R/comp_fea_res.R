##' Dot plot for comparing the top ranking functional categories from different 
##' functional enrichment analysis (FEA) results. The functional categories are 
##' plotted in the order defined by their mean rank across the corresponding 
##' FEA results.
##' 
##' The `comp_fea_res` function computes the mean rank for each functional 
##' category across different FEA result instances and then re-ranks them based 
##' on that.Since the functional categories are not always present in all 
##' enrichment results, the mean rank of a functional category is corrected by 
##' an adjustment factor that is the number of enrichment result methods used 
##' divided by the number of occurences of a functional category. For instance, 
##' if a functional category is only present in the result of one method, its 
##' mean rank will be increased accordingly. Subsequently, the re-ranked 
##' functional categories are compared in a dot plot where the colors represent 
##' the values of the enrichment statistic chosen under the \code{rank_stat} 
##' argument.
##'
##' @title Plot for Comparing Ranking Results of FEA Methods
##' @param table_list a named list of tibbles extracted from feaResult objects,
##' e.g. generated with different FEA methods.
##' @param rank_stat character(1), column name of the enrichment statisic used 
##' for ranking the functional categories, e.g. 'pvalue' or 'p.adjust'. Note,
##' the chosen column name needs to be present in each tibble of 'table_list'.
##' @param Nshow integer defining the number of the top functional categories 
##' to display in the plot after re-ranking them across FEA methods
##' @importFrom ggplot2 aes_
##' @importFrom DOSE theme_dose
##' @return ggplot2 graphics object
##' @examples 
##' data(drugs)
##' \dontrun{
##' table_list <- list("dup_hyperG" = result(dup_hyperG_res), 
##'                   "mGSEA" = result(mgsea_res), 
##'                   "mabs" = result(mabs_res), 
##'                   "hyperG" = result(hyperG_res), 
##'                   "GSEA" = result(gsea_res))
##' comp_fea_res(table_list, rank_stat="pvalue", Nshow=20)
##' }
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
    cat_mrk_adj <- vapply(cat_ranks_list, function(i){
        m <- mean(i)*5/length(i)
    }, FUN.VALUE = numeric(1))
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


