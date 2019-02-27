cor_sig_search <- function(query, refdb, method){
  res_list <- NULL
  # make sure rownames of query and refdb are exactly same
  common_gene <- intersect(rownames(query), rownames(refdb))
  query2 <- as.matrix(query[common_gene,])
  colnames(query2) <- colnames(query)
  refdb <- refdb[common_gene,]
  for(i in 1:ncol(query2)){
    cor <- as.numeric(cor(query2[,i], refdb, method = method))
    names(cor) <- colnames(refdb)
    trend=cor
    trend[cor>=0]="up"
    trend[cor<0]="down"
    res <- data.frame(set=names(cor), trend=trend, cor_score = cor, stringsAsFactors = FALSE)
    res_list <- c(res_list, list(res))
  }
  names(res_list) <- colnames(query)
  if(length(res_list)==1) return(res_list[[1]])
  return(res_list)
}

#' Gene expression signature search with correlation coefficient
#' @title gess_cor
#' @param qSig `qSig` object, The 'gess_method' slot should be 'Cor'. The ‘qsig’ slot could be a matrix representing gene expression values, 
#' the reference database could also store genome-wide expression values of treatment samples
#' @param method a character string indicating which correlation coefficient (or covariance) is to be computed. 
#' One of "spearman" (default), "kendall", or "pearson": can be abbreviated.
#' @param chunk_size size of chunk per processing
#' @return gessResult object
#' @export
gess_cor <- function(qSig, method, chunk_size=5000){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  #stopifnot(validObject(qSig))
  if(qSig@gess_method != "Cor"){
    stop("The 'gess_method' slot of 'qSig' should be 'Cor' if using 'gess_cor' function")
  }
  query <- qSig@qsig
  se <- qSig@refdb
  dmat <- assay(se)
  ceil <- ceiling(ncol(dmat)/chunk_size)
  resultDF <- data.frame()
  for(i in 1:ceil){
    dmat_sub <- dmat[,(chunk_size*(i-1)+1):min(chunk_size*i, ncol(dmat))]
    mat <- as(dmat_sub, "matrix")
    cor_res <- cor_sig_search(query=query, refdb=mat, method=method)
    resultDF <- rbind(resultDF, data.frame(cor_res))
  }
  set = resultDF$set
  new <- as.data.frame(t(sapply(1:length(set), function(i)
    unlist(strsplit(as.character(set[i]), "__")))), stringsAsFactors=FALSE)
  colnames(new) = c("pert", "cell", "type")
  resultDF <- data.frame(new, resultDF[,-1])
  resultDF <- resultDF[order(abs(resultDF$cor_score), decreasing = TRUE), ]
  row.names(resultDF) <- NULL
  # add target column
  target <- suppressMessages(get_targets(resultDF$pert))
  res <- left_join(resultDF, target, by=c("pert"="drug_name"))
  
  x <- gessResult(result = as_tibble(res),
                  qsig = qSig@qsig,
                  gess_method = qSig@gess_method,
                  refdb = qSig@refdb)
  return(x)
}
