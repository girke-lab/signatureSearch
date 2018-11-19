#' Summarize ranks of drugs in GESS result across cell types
#' @title cell_rank_sum
#' @param gessResult `gessResult` object
#' @return tibble
#' @importFrom dplyr bind_cols
#' @export
cell_rank_sum <- function(gessResult){
  if(!is(gessResult, "gessResult")) stop("The 'gessResult' should be an object of 'gessResult' class")
  tb <- result(gessResult)[,c(1:2)]
  tb <- bind_cols(rank = 1:nrow(tb), tb)
  # rank_sum <- matrix(NA, nrow = length(unique(tb$pert)), ncol = length(unique(tb$cell)))
  # rownames(rank_sum)=as.character(unique(tb$pert))
  # colnames(rank_sum)=as.character(unique(tb$cell))
  # for(i in 1:nrow(rank_sum)){
  #   for(j in 1:ncol(rank_sum)){
  #     rank_sum[i,j]=as.numeric(tb[tb$pert==rownames(rank_sum)[i]&tb$cell==colnames(rank_sum)[j], "rank"])
  #   }
  # }
  dl <- split(as.data.frame(tb)[,c("rank","cell")], tb$pert)
  cell_name=unique(tb$cell)
  dl_num <- sapply(dl, function(x){
    num <- x$rank
    names(num) <- as.character(x$cell)
    num2 <- num[cell_name]
    names(num2) <- cell_name
    return(num2)
  }, simplify = FALSE)
  mat <- t(as.data.frame(dl_num, check.names=FALSE))
  mat <- as.data.frame(mat)
  min <- apply(mat,1,min, na.rm=TRUE)
  max <- apply(mat,1,max, na.rm=TRUE)
  res <- data.frame(pert=rownames(mat), mat, min=min, max=max)
  res <- res[order(res$min),]
  return(res)
}

#' Append two columns (score_column_grp1, score_column_grp2) to GESS result
#' @title sim_score_grp
#' @param tib tibble object from result(gessResult)
#' @param grp1 character vector, group 1 of cell types, e.g., tumor cell types
#' @param grp2 character vector, group 2 of cell types, e.g., normal cell types
#' @param score_column character, column name of similairity scores to be grouped 
#' @return tibble
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
