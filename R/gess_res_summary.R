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
