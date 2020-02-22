##' Dot plot for comparing the top ranking functional categories from different 
##' functional enrichment analysis (FEA) results. The functional categories are 
##' plotted in the order defined by their mean rank across the corresponding 
##' FEA results.
##' 
##' The `comp_fea_res` function computes the mean rank for each functional 
##' category across different FEA result instances and then re-ranks them based 
##' on that. Since the functional categories are not always present in all 
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
##' @param Nchar integer defining number of characters displayed (exceeded 
##' characters were replaced by '...') in the description of each item
##' @param scien TRUE or FALSE, indicating whether the rank_stat is rounded to
##' the scientific format with 3 digits
##' @param ... Other arguments passed on to \code{\link[ggplot2]{geom_point}}
##' @importFrom ggplot2 aes_
##' @importFrom DOSE theme_dose
##' @return ggplot2 graphics object
##' @examples 
##' method1 <- data.frame("ID"=paste0("GO:", 1:5), 
##'                       "Description"=paste0("desc", 1:5),
##'                       "pvalue"=c(0.0001, 0.002, 0.004, 0.01, 0.05))
##' method2 <- data.frame("ID"=paste0("GO:", c(1,3,5,4,6)), 
##'                       "Description"=paste0("desc", c(1,3,5,4,6)),
##'                       "pvalue"=c(0.0003, 0.0007, 0.003, 0.006, 0.04))
##' table_list <- list("method1" = method1, "method2"=method2) 
##' comp_fea_res(table_list, rank_stat="pvalue", Nshow=20)
##' @export

comp_fea_res <- function(table_list, rank_stat="pvalue", Nshow=20, 
                         Nchar=50, scien=FALSE, ...){
    if(is.null(names(table_list))){
        stop(paste('The "table_list" should be a list with names,', 
             'which could be names of the enrichment methods'))
    }
    newtb <- function(i){
        tb <- table_list[[i]]
        tb_order <- tb[order(tb[[rank_stat]]),]
        tb2 <- data.frame(tb_order[, c("Description", rank_stat)], 
                          Method=names(table_list)[i], 
                          rank=seq_len(nrow(tb_order)))
        return(tb2)
    }
    df <- do.call(rbind, lapply(seq_along(table_list), newtb))
    df$Description <- vec_char_redu(df$Description, Nchar=Nchar)
    if(scien){
        df[[rank_stat]] <- signif(df[[rank_stat]], 3)
    }
    # re-ranking according to mean ranks and adjusted for number of support
    # methods (*length(table_list)/number of support methods)
    cat_ranks_list <- split(df$rank, df$Description) 
    cat_mrk_adj <- vapply(cat_ranks_list, function(i){
        m <- mean(i)*length(table_list)/length(i)
    }, FUN.VALUE = numeric(1))
    cat_rerk <- cat_mrk_adj[order(cat_mrk_adj)]
    cat_top <- names(cat_rerk)[seq_len(Nshow)]
    df2 <- df[df$Description %in% cat_top,]
    df2$Description <- factor(df2$Description, levels = rev(cat_top), 
                            ordered = TRUE)
    p <- ggplot(df2, aes_(x = ~Method, y = ~Description)) +
    geom_point(size=4, ...) + aes_string(color=rank_stat) + 
    scale_colour_gradient(low="red", high="blue") +
    theme_dose(font.size=12) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
    return(p)
}


