sp_sig_search <- function(query, refdb){
  res_list <- NULL
  # make sure rownames of query and refdb are exactly same
  common_gene <- intersect(rownames(query), rownames(refdb))
  query2 <- as.matrix(query[common_gene,])
  colnames(query2) <- colnames(query)
  refdb <- refdb[common_gene,]
  for(i in 1:ncol(query2)){
    cor <- cor(query2[,i], refdb, method = "spearman")
    names(cor) <- colnames(refdb)
    trend=cor
    trend[cor>=0]="up"
    trend[cor<0]="down"
    res <- data.frame(set=names(cor), trend=trend, sp_score = cor, stringsAsFactors = FALSE)
    res_list <- c(res_list, list(res))
  }
  names(res_list) <- colnames(query)
  if(length(res_list)==1) return(res_list[[1]])
  return(res_list)
}

#' Gene expression signature search with Spearman correlation
#' @title gess_sp
#' @param qSig `qSig` object, The 'gess_method' slot of 'qSig' should be 'SP'
#' @param chunk_size size of chunk per processing
#' @return gessResult object
#' @export
gess_sp <- function(qSig, chunk_size=5000){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  #stopifnot(validObject(qSig))
  if(qSig@gess_method != "SP"){
    stop("The 'gess_method' slot of 'qSig' should be 'SP' if using 'gess_sp' function")
  }
  query <- qSig@qsig
  se <- qSig@refdb
  dmat <- assay(se)
  ceil <- ceiling(ncol(dmat)/chunk_size)
  resultDF <- data.frame()
  for(i in 1:ceil){
    dmat_sub <- dmat[,(chunk_size*(i-1)+1):min(chunk_size*i, ncol(dmat))]
    mat <- as(dmat_sub, "matrix")
    sp_res <- sp_sig_search(query=query, refdb=mat)
    resultDF <- rbind(resultDF, data.frame(sp_res))
  }
  resultDF <- resultDF[order(abs(resultDF$sp_score), decreasing = TRUE), ]
  row.names(resultDF) <- NULL
  x <- gessResult(result = as_tibble(resultDF),
                  qsig = qSig@qsig,
                  gess_method = qSig@gess_method,
                  refdb = qSig@refdb)
  return(x)
}
