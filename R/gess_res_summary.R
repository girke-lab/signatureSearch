#' For drugs in GESS result, get their ranks in different cell types
#' @title Summarize ranks of drugs in different cells
#' @param gessResult `gessResult` object
#' @return data.frame
#' @importFrom dplyr bind_cols
#' @examples 
#' db_path <- system.file("extdata", "sample_db.h5", 
#'                        package = "signatureSearch")
#' library(signatureSearchData)
#' sample_db <- readHDF5chunk(db_path, colindex=1:100)
#' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
#' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
#' query = as.numeric(query_mat); names(query) = rownames(query_mat)
#' upset <- head(names(query[order(-query)]), 150)
#' downset <- tail(names(query[order(-query)]), 150)
#' qsig_cmap <- qSig(query = list(upset=upset, downset=downset), 
#'                   gess_method = "CMAP", refdb = db_path)
#' cmap <- gess_cmap(qSig=qsig_cmap, chunk_size=5000)
#' df <- drug_cell_ranks(cmap)
#' @export
drug_cell_ranks <- function(gessResult){
  if(!is(gessResult, "gessResult")) 
    stop("The 'gessResult' should be an object of 'gessResult' class")
  tb <- result(gessResult)[,c(seq_len(2))]
  tb <- bind_cols(rank = seq_len(nrow(tb)), tb)
  dl <- split(as.data.frame(tb)[,c("rank","cell")], tb$pert)
  cell_name=unique(tb$cell)
  dl_num <- lapply(dl, function(x){
    num <- x$rank
    names(num) <- as.character(x$cell)
    num2 <- num[cell_name]
    names(num2) <- cell_name
    return(num2)})
  mat <- t(as.data.frame(dl_num, check.names=FALSE))
  mat <- as.data.frame(mat)
  min <- apply(mat,1,min, na.rm=TRUE)
  mean <- apply(mat,1,mean, na.rm=TRUE)
  max <- apply(mat,1,max, na.rm=TRUE)
  res <- data.frame(pert=rownames(mat), mat, min=min, mean=mean, max=max)
  res <- res[order(res$min),]
  rownames(res)=NULL
  return(res)
}

#' Append two columns (score_column_grp1, score_column_grp2) to GESS result,
#' which represent cell summary scores of drugs in a cell group.
#' 
#' @title Get summary scores of drugs in cell groups
#' @param tib tibble in gessResult object. 
#' @param grp1 character vector, group 1 of cell types, e.g., tumor cell types
#' @param grp2 character vector, group 2 of cell types, e.g., normal cell types
#' @param score_column character, column name of similarity scores to be 
#' grouped 
#' @return tibble
#' @examples 
#' db_path <- system.file("extdata", "sample_db.h5", 
#'                        package = "signatureSearch")
#' library(signatureSearchData)
#' sample_db <- readHDF5chunk(db_path, colindex=1:100)
#' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
#' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
#' query = as.numeric(query_mat); names(query) = rownames(query_mat)
#' upset <- head(names(query[order(-query)]), 150)
#' downset <- tail(names(query[order(-query)]), 150)
#' qsig_lincs <- qSig(query = list(upset=upset, downset=downset), 
#'                    gess_method = "LINCS", refdb = db_path)
#' lincs <- gess_lincs(qsig_lincs, sortby="NCS", tau=FALSE)
#' df <- sim_score_grp(result(lincs), grp1="SKB", grp2="MCF7", "NCS")
#' @export

sim_score_grp <- function(tib, grp1, grp2, score_column){
  ## Summary across group of cell lines
  ctgrouping <- paste(tib$pert, tib$type, sep="__")
  
  tib_grp1 <- dplyr::filter(tib, cell %in% grp1)
  ctgrouping1 <- paste(tib_grp1$pert, tib_grp1$type, sep="__")
  cs1 <- tib_grp1[[score_column]]
  qmax1 <- tapply(cs1, ctgrouping1, function(x) { 
    q <- quantile(x, probs=c(0.33, 0.67))
    ifelse(abs(q[2]) >= abs(q[1]), q[2], q[1])
  })
  qmax1 <- qmax1[ctgrouping]
  
  tib_grp2 <- filter(tib, cell %in% grp2)
  ctgrouping2 <- paste(tib_grp2$pert, tib_grp2$type, sep="__")
  cs2 <- tib_grp2[[score_column]]
  qmax2 <- tapply(cs2, ctgrouping2, function(x) { 
    q <- quantile(x, probs=c(0.33, 0.67))
    ifelse(abs(q[2]) >= abs(q[1]), q[2], q[1])
  })
  qmax2 <- qmax2[ctgrouping]
  name_grp1 <- paste0(score_column, "_grp1")
  name_grp2 <- paste0(score_column, "_grp2")
  cname <- colnames(tib)
  tib %<>% dplyr::mutate(qmax1, qmax2)
  colnames(tib) <- c(cname, name_grp1, name_grp2)
  return(tib)
}

#' The GESS result from a comprehensive signature database, e.g., LINCS,
#' contains similarity scores of drugs in different cell types, 
#' including normal and tumor. The similarity scores of a drug in cell types 
#' can be summarized and visualized. The summary method is according to the 
#' "Summarization Across Cell Lines" method from Subramanian et al., 2017.
#'
#' @title GESS result visualization
#' @param gess_tb tibble in \code{\link{gessResult}} object, can be accessed
#' via \code{\link{result}} method.
#' @param drugs character vector, a list of interesting drugs
#' @param col name of the score column in 'gess_tb', e.g., "NCS" 
#' @return plot
#' @references  
#' Subramanian, A., Narayan, R., Corsello, S. M., Peck, D. D., 
#' Natoli, T. E., Lu, X., … Golub, T. R. (2017). A Next Generation 
#' Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell, 
#' 171(6), 1437–1452.e17. \url{https://doi.org/10.1016/j.cell.2017.10.049}
#' @importFrom readr read_tsv
#' @importFrom dplyr mutate
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @examples 
#' db_path <- system.file("extdata", "sample_db.h5", 
#'                        package = "signatureSearch")
#' library(signatureSearchData)
#' sample_db <- readHDF5chunk(db_path, colindex=1:100)
#' ## get "vorinostat__SKB__trt_cp" signature drawn from sample databass
#' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
#' query = as.numeric(query_mat); names(query) = rownames(query_mat)
#' upset <- head(names(query[order(-query)]), 150)
#' downset <- tail(names(query[order(-query)]), 150)
#' qsig_lincs <- qSig(query = list(upset=upset, downset=downset), 
#'                    gess_method = "LINCS", refdb = db_path)
#' lincs <- gess_lincs(qsig_lincs, sortby="NCS", tau=FALSE)
#' data(drugs)
#' gess_res_vis(result(lincs), drugs, col="NCS")
#' @export
gess_res_vis <- function(gess_tb, drugs, col){
  ext_path  <- system.file("extdata", package="signatureSearch")
  cell_path <- paste0(ext_path,"/cell_info.tsv")
  cell_info <- suppressMessages(read_tsv(cell_path))
  cells = unique(gess_tb$cell)
  
  # Summarize NCS across groups
  tumor <- cell_info %>% filter(cell_type=="tumor")
  tumor = intersect(as.character(tumor$cell_id), cells)
  normal <- cell_info %>% filter(cell_type=="normal")
  normal = intersect(as.character(normal$cell_id), cells)
  tb_grp <- sim_score_grp(gess_tb, tumor, normal, score_column=col)

  data1 = gess_tb %>% mutate(rank=seq_len(nrow(gess_tb))) %>% 
    filter(pert %in% drugs) %>% mutate(pert = factor(pert, levels=drugs)) %>% 
    left_join(cell_info[,c("cell_id","cell_type")], by=c("cell"="cell_id"))
  
  data2 = tb_grp %>% filter(pert %in% drugs) %>% 
    mutate(pert = factor(pert, levels=drugs))
  data2 = data2[,c("pert", paste0(col, "_grp1"), paste0(col,"_grp2"))]
  colnames(data2) = c("pert", "tumor", "normal")
  data2 %<>% distinct() %>%
    reshape2::melt(id.vars=c("pert"), variable.name = "cell_type", 
                   value.name=paste0(col, "_grp"))
  
  p = ggplot(data1) + 
    geom_point(data=data1, 
               aes_string(x="pert", y=col, fill = "cell", shape = "cell_type", 
                          colour = "cell"), size=2.5) +
    geom_point(data=data2, 
               aes_string("pert", paste0(col, "_grp"), shape = "cell_type"), 
               size=2.5) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank())
  p
}
cell_type = cell_id = pert = NULL
