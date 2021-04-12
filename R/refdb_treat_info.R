#' Get treatment information including perturbation name, cell type and 
#' perturbation type from the reference database
#' 
#' @title Get Treatment Information
#' @param refdb character(1), one of "lincs", "lincs_expr", "cmap" or "cmap_expr"
#' when using the pre-generated CMAP/LINCS databases or path to the HDF5 file
#' generated with the \code{\link{build_custom_db}} function. The details is 
#' shown in the 'refdb' argument of the \code{\link{qSig}} function
#' @return tibble object with 'pert', 'cell', 'pert_type' columns 
#' @examples 
#' refdb <- system.file("extdata", "sample_db.h5", package="signatureSearch")
#' treat_info <- get_treat_info(refdb)
#' @export
get_treat_info <- function(refdb){
    db_path <- determine_refdb(refdb)
    treats <- as.character(HDF5Array(db_path, name="colnames"))
    treat_df <- as.data.frame(t(vapply(seq_along(treats), function(i)
        unlist(strsplit(as.character(treats[i]), "__")), 
        FUN.VALUE=character(3))), stringsAsFactors=FALSE)
    colnames(treat_df) <- c("pert", "cell", "pert_type")
    treat_df <- as_tibble(treat_df)
    return(treat_df)
}

#' Bar plot of number of perturbations/compounds tested in cell types where 
#' cell types are grouped by 'primary site'.
#' 
#' @title Number of Tests in Cell Types
#' @inheritParams get_treat_info
#' @return Faceted bar plot 
#' @examples
#' refdb <- system.file("extdata", "sample_db.h5", package="signatureSearch")
#' cellNtestPlot(refdb)
#' @export
cellNtestPlot <- function(refdb){
    data("cell_info", envir=environment())
    cell_site <- cell_info$primary_site
    names(cell_site) <- cell_info$cell_id
    
    treat_df <- get_treat_info(refdb)
    ntest <- table(treat_df$cell)
    df <- data.frame(Ntest=as.numeric(ntest), Cell=names(ntest), 
                     PrimarySite=cell_site[names(ntest)])
    
    p <- ggplot(df, aes(x=Cell, y=Ntest, fill=Cell)) + 
        geom_bar(stat="identity", width=1) +
        geom_text(aes(label=Ntest), vjust=0.5, hjust=0.05, size=3) +
        coord_flip() +
        facet_grid(PrimarySite~., switch="y", scales="free", space="free") +
        theme(panel.spacing=unit(0.1, "lines"), legend.position="none",
              strip.text.y = element_text(size=6))
    p
}

## get rid of "Undefined global functions or variables" note
globalVariables(c("Cell", "Ntest")) 
