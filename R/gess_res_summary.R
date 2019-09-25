#' The \code{drug_cell_ranks} function returns from a \code{gessResult} object 
#' the ranks of the perturbagens (e.g. drugs) for each cell type. The results 
#' are arranged in separate columns of a \code{data.frame}. Additionally, it 
#' includes in the last columns summary ranking statistics across all cell 
#' types, such as min, mean and max values.
#' @title Summary ranking statistics across cell types
#' @param gessResult `gessResult` object
#' @return data.frame
#' @importFrom dplyr bind_cols
#' @examples 
#' gr <- gessResult(result=dplyr::tibble(pert=c("p1", "p1", "p2", "p3"),
#'                                       cell=c("MCF7", "SKB", "MCF7", "SKB"),
#'                                       type=rep("trt_cp", 4),
#'                                       NCS=c(1.2, 1, 0.9, 0.6)),
#'                  query=list(up="a", down="b"), 
#'                  gess_method="LINCS", refdb="path/to/refdb")
#' df <- drug_cell_ranks(gr)
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

#' Function appends two columns (score_column_grp1, score_column_grp2) to GESS 
#' result tibble. The appended columns contain summary scores for groups of 
#' cell types, such as normal and tumor cells.
#' @title Summary Scores by Groups of Cell Types
#' @param tib tibble in gessResult object 
#' @param grp1 character vector, group 1 of cell types, e.g., tumor cell types
#' @param grp2 character vector, group 2 of cell types, e.g., normal cell types
#' @param score_column character, column name of similarity scores to be 
#' grouped 
#' @return tibble
#' @examples 
#' gr <- gessResult(result=dplyr::tibble(pert=c("p1", "p1", "p2", "p3"),
#'                                       cell=c("MCF7", "SKB", "MCF7", "SKB"),
#'                                       type=rep("trt_cp", 4),
#'                                       NCS=c(1.2, 1, 0.9, 0.6)),
#'                  query=list(up="a", down="b"), 
#'                  gess_method="LINCS", refdb="path/to/refdb")
#' df <- sim_score_grp(result(gr), grp1="SKB", grp2="MCF7", "NCS")
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

#' The function allows to summarize the ranking scores of selected
#' perturbagens for GESS results across cell types along with cell type 
#' classifications, such as normal and tumor cells. In the resulting plot
#' the perturbagens are drugs (along x-axis) and the ranking scores are LINCS'
#' NCS values (y-axis). For each drug the NCS values are plotted for each cell
#' type as differently colored dots, while their shape indicates the cell type
#' class.
#'
#' @title GESS Result Visualization
#' @param gess_tb tibble in the 'result' slot of the \code{\link{gessResult}} 
#' object, can be extracted via \code{\link{result}} accessor function
#' @param drugs character vector of selected drugs
#' @param col character(1), name of the score column in 'gess_tb', e.g., "NCS"
#' if the result table is from LINCS method. Can also be set as "rank", 
#' this way it will show the ranks of each drug in different cell types.
#' @param cell_group character(1), one of "all", "normal", or "tumor". 
#' If "all", it will show scores of each drug in both tumor and normal cell 
#' types. If "normal" or "tumor", it will only show normal or tumor cell types.
#' @return plot visualizing GESS results
#' @references  
#' Subramanian, A., Narayan, R., Corsello, S. M., Peck, D. D., 
#' Natoli, T. E., Lu, X., Golub, T. R. (2017). A Next Generation 
#' Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell, 
#' 171 (6), 1437-1452.e17. URL: https://doi.org/10.1016/j.cell.2017.10.049
#' @importFrom readr read_tsv
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_colour_gradient
#' @importFrom utils data
#' @examples 
#' gr <- gessResult(result=dplyr::tibble(pert=c("p1", "p1", "p2", "p3"),
#'                                       cell=c("MCF7", "SKB", "MCF7", "SKB"),
#'                                       type=rep("trt_cp", 4),
#'                                       NCS=c(1.2, 1, 0.9, 0.6)),
#'                  query=list(up="a", down="b"), 
#'                  gess_method="LINCS", refdb="path/to/refdb")
#' gess_res_vis(result(gr), drugs=c("p1","p2"), col="NCS")
#' @export
gess_res_vis <- function(gess_tb, drugs, col, cell_group="all"){
  # ext_path  <- system.file("extdata", package="signatureSearch")
  # cell_path <- paste0(ext_path,"/cell_info.tsv")
  # cell_info <- suppressMessages(read_tsv(cell_path))
  data("cell_info", envir=environment())
  cells = unique(gess_tb$cell)
  if(col=="rank"){
      gess_tb <- data.frame(gess_tb, rank=seq_len(dim(gess_tb)[1]))
  }
  # Summarize NCS across groups
  tumor <- cell_info %>% filter(cell_type=="tumor")
  tumor = intersect(as.character(tumor$cell_id), cells)
  normal <- cell_info %>% filter(cell_type=="normal")
  normal = intersect(as.character(normal$cell_id), cells)
  tb_grp <- sim_score_grp(gess_tb, tumor, normal, score_column=col)

  data1 = gess_tb %>% 
    filter(pert %in% drugs) %>% mutate(pert = factor(pert, levels=drugs)) %>% 
    left_join(cell_info[,c("cell_id","cell_type")], by=c("cell"="cell_id")) %>%
    rename("cell_class"= cell_type) %>%
    rename("cell_type"=cell)
  
  data2 = tb_grp %>% filter(pert %in% drugs) %>% 
    mutate(pert = factor(pert, levels=drugs))
  data2 = data2[,c("pert", paste0(col, "_grp1"), paste0(col,"_grp2"))]
  colnames(data2) = c("pert", "tumor", "normal")
  data2 %<>% distinct() %>%
    reshape2::melt(id.vars=c("pert"), variable.name = "cell_class", 
                   value.name=paste0(col, "_grp"))
  if(cell_group != "all"){
      data1 %<>% filter(cell_type == cell_group)
      data2 %<>% filter(cell_type == cell_group) 
  }
  p = ggplot(data1) + 
    geom_point(data=data1, 
               aes_string(x="pert", y=col, fill="cell_type", shape="cell_class",
                          colour = "cell_type"), size=2.5) +
    geom_point(data=na.omit(data2), 
               aes_string("pert", paste0(col, "_grp"), shape = "cell_class"), 
               size=2.5) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank())
  p
}

globalVariables(c("cell_type", "cell_id", "pert", "cell_info")) 
