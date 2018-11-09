#' @importFrom utils read.delim
#' @importFrom stats quantile
.lincsScores <- function(esout, upset, downset, minTauRefSize, ES_NULL="Default", taurefList="Default") {
  ## P-value and FDR for WTCS based on ESnull from random queries where p-value = sum(ESrand > ES_obs)/Nrand
  if(ES_NULL != "Default") {WTCSnull <- read.delim(ES_NULL)} else {
    ext_path <- system.file("extdata", package="signatureSearch")
    if(! file.exists(file.path(ext_path, "ES_NULL.txt"))){
      download.file("http://biocluster.ucr.edu/~yduan004/LINCS_db/ES_NULL.txt", file.path(ext_path, "ES_NULL.txt"), quiet = TRUE)
    }
    WTCSnull <- read.delim(file.path(ext_path, "ES_NULL.txt")) 
  }
  WTCSnull[WTCSnull[, "Freq"]==0, "Freq"] <- 1 # Add pseudo count of 1 where Freq is zero 
  myrounding <- max(nchar(as.character(WTCSnull[,"WTCS"]))) - 3 # Three because of dot and minus sign
  es_round <- round(as.numeric(esout), myrounding) # Assures same rounding used for WTCSnull computation
  WTCS_pval <- sapply(es_round, function(x) {
    sum(WTCSnull[abs(WTCSnull[,"WTCS"]) > abs(x), "Freq"]) / sum(WTCSnull[,"Freq"])	
  })
  WTCS_fdr <- p.adjust(WTCS_pval, "fdr")
  ## Normalized connectivity score (NCS)
  grouping <- paste(gsub("^.*?__", "", names(esout)), as.character(ifelse(esout > 0, "up", "down")), sep="__")
  es_na <- as.numeric(esout)
  es_na[es_na == 0] <- NA # eliminates zeros from mean calculation; zeros have high impact on NCS values due to their high frequency
  groupmean <- tapply(es_na, grouping, mean, na.rm=TRUE)
  groupmean[is.na(groupmean)] <- 0 # In case groups contain only zeros, NA/NaN values are introduced in mean calculation which are reset to zeros here
  groupmean[groupmean==0] <- 10^-12 # Set zeros (can be from non NAs) to small value to avoid devision by zero 
  ncs <- as.numeric(esout) / abs(groupmean[grouping]) # without abs() sign of neg values would switch to pos
  ## Tau calculation requires reference NCS lookup DB
  ## performs: sign(ncs_query) * 100/N sum(abs(ncs_ref) < abs(ncs_query))
  if(! file.exists(file.path(ext_path, "taurefList.rds"))){
    download.file("http://biocluster.ucr.edu/~yduan004/LINCS_db/taurefList.rds", file.path(ext_path, "taurefList.rds"), quiet = TRUE)
  }
  taurefList9264 <- readRDS(file.path(ext_path, "taurefList.rds"))
  ncs_query <- ncs; names(ncs_query) <- names(esout)
  queryDB_refDB_match <- unique(unlist(lapply(taurefList9264, rownames))) %in% names(ncs_query)
  if(!all(queryDB_refDB_match)) warning(paste0("QueryDB and tauRefDB differ by ", round(100 * sum(!queryDB_refDB_match)/length(queryDB_refDB_match),1), "% of their entries. Accurate tau computation requires close to 0% divergence."))
  ncs_query_list <- split(ncs_query, factor(gsub("^.*?__", "", names(ncs_query))))
  tau <- lapply(names(ncs_query_list), function(x) {
    tmpDF <- taurefList9264[[x]]
    ncs_query_match <- names(ncs_query_list[[x]])[names(ncs_query_list[[x]]) %in% rownames(tmpDF)]
    if(length(ncs_query_match)>0) {
      tmpDF <- tmpDF[ncs_query_match, , drop=FALSE]	
      # sign(ncs_query_list[[x]]) * 100/ncol(tmpDF) * rowSums(abs(tmpDF) < abs(ncs_query_list[[x]]))
      sign(ncs_query_list[[x]]) * 100/ncol(tmpDF) * rowSums(abs(tmpDF) < abs(round(ncs_query_list[[x]], 2))) # rounded as in ref db
    } else {
      NULL
    }
  })
  tau <- unlist(tau)
  tau <- tau[names(ncs_query)]
  tauRefSize <- sapply(taurefList9264, ncol)[gsub("^.*?__", "", names(tau))]
  tau[tauRefSize < minTauRefSize] <- NA
  
  ## Summary across cell lines (NCSct)
  ctgrouping <- gsub("__.*__", "__", names(esout))
  qmax <- tapply(ncs, ctgrouping, function(x) { 
    q <- quantile(x, probs=c(0.33, 0.67))
    ifelse(abs(q[2]) >= abs(q[1]), q[2], q[1])
  })
  qmax <- qmax[ctgrouping]
  ## Organize result in data.frame
  new <- as.data.frame(t(sapply(1:length(esout), function(i)
    unlist(strsplit(as.character(names(esout)[i]), "__")))), stringsAsFactors=FALSE)
  colnames(new) = c("pert", "cell", "type")
  resultDF <- data.frame(new, 
                         trend = as.character(ifelse(esout > 0, "up", "down")),
                         WTCS = as.numeric(esout),
                         WTCS_Pval = WTCS_pval,
                         WTCS_FDR = WTCS_fdr,
                         NCS = ncs,
                         Tau = tau,
                         TauRefSize=tauRefSize,
                         NCSct = qmax,
                         NCSct_group=ctgrouping,
                         N_upset = length(upset),
                         N_downset = length(downset), stringsAsFactors = FALSE)
  row.names(resultDF) <- NULL
  return(resultDF)
}


## Define enrichment function according to Subramanian et al, 2005
## Note: query corresponds to gene set, here Q.
.enrichScore <- function(sigvec, Q, type) {
  ## Preprocess arguments
  L <- names(sigvec)
  R <- as.numeric(sigvec)
  N <- length(L)
  NH <- length(Q)
  Ns <- N - NH
  hit_index <- as.numeric(L %in% Q)	
  miss_index <- 1 - hit_index	
  R <- abs(R^type)
  ## Compute ES
  NR <- sum(R[hit_index == 1])
  ESvec <- cumsum((hit_index * R * 1/NR) - (miss_index * 1/Ns)) 
  ES <- ESvec[which.max(abs(ESvec))]
  return(ES)
}

#' lincsEnrich
#' @param se SummarizedExperiment object
#' @param upset character vector represents up regulated gene sets
#' @param downset character vector represents down regulated gene sets
#' @param sortby by which the result data frame is sorted
#' @param type exponent of GSEA algorithm
#' @param output one of `esonly` and `all`.
#' @param ES_NULL path to the ES_NULL file. ES_null distribution is generated with random queryies for computing nominal P-values for ES by 
#' `randQueryES_slurm` or `randQueryES_local` function. If `ES_NULL` is set as `Default`, it uses the ES_null distribution that we generated.
#' @param taurefList path to the "taurefList.rds" file. If set as `Default`, it uses the "taurefList.rds" file that we generated.
#' @param minTauRefSize minimum size of reference data to compute Tau score
#' @param chunk_size size of chunk per processing
#' @return data.frame
#' @importFrom DelayedArray apply
#' @export
lincsEnrich <- function(se, upset, downset, sortby="NCS", type=1, output="all", ES_NULL="Default", taurefList="Default", minTauRefSize=500, chunk_size=5000) {
  mycolnames <- c("WTCS", "NCS", "Tau", "NCSct", "N_upset", "N_downset", NA)
  if(!any(mycolnames %in% sortby)) stop("Unsupported value assinged to sortby.")
  
  ## Run .enrichScore on upset and downset
  ## When both upset and downset are provided 
  dmat <- assay(se)
  
  ceil <- ceiling(ncol(dmat)/chunk_size)
  ESout=NULL
  for(i in 1:ceil){
    dmat_sub <- dmat[,(chunk_size*(i-1)+1):min(chunk_size*i, ncol(dmat))]
    mat <- as(dmat_sub, "matrix")
    if(length(upset)>0 & length(downset)>0) {
      ESup <- apply(mat, 2, function(x) .enrichScore(sigvec=sort(x, decreasing = TRUE), Q=upset, type=type))
      ESdown <- apply(mat, 2, function(x) .enrichScore(sigvec=sort(x, decreasing = TRUE), Q=downset, type=type)) 
      ESout1 <- ifelse(sign(ESup) != sign(ESdown), (ESup - ESdown)/2, 0)
      ## When only upset is provided
    } else if(length(upset)>0 & length(downset)==0) {
      ESup <- apply(mat, 2, function(x) .enrichScore(sigvec=sort(x, decreasing = TRUE), Q=upset, type=type))
      ESout1 <- ESup
      ## When only downset is provided
    } else if(length(upset)==0 & length(downset)>0) {
      ESdown <- apply(mat, 2, function(x) .enrichScore(sigvec=sort(x, decreasing = TRUE), Q=downset, type=type))
      ESout1 <- -ESdown
      ## When none are provided (excluded by input validity check already!)
    }
    ESout <- c(ESout, ESout1)
  }
  
  ## Assmble output 
  if(output=="esonly") {	
    return(ESout)
  }
  if(output=="all") {
    resultDF <- .lincsScores(esout=ESout, upset=upset, downset=downset, minTauRefSize=minTauRefSize, ES_NULL=ES_NULL, taurefList=taurefList)
  }
  if(!is.na(sortby)) {
    resultDF <- resultDF[order(abs(resultDF[,sortby]), decreasing=TRUE), ]
  } else {
    resultDF <- resultDF
  }
  row.names(resultDF) <- NULL
  return(resultDF)
}

#' LINCS method for GESS
#' 
#' @title gess_lincs
#' @param qSig `qSig` object, The 'gess_method' slot should be 'LINCS'
#' @param ES_NULL path to the ES_NULL file. ES null distribution is generated with random queries for computing nominal P-values for ES by 
#' `randQueryES_slurm` or `randQueryES_local` function. 
#' If `ES_NULL` is set as `Default`, it uses the ES_null distribution that we generated against LINCS database.
#' @param taurefList path to the "taurefList.rds" file generated with `queryReferenceDB` function. 
#' If set as `Default`, it uses the "taurefList.rds" file that we generated against LINCS database.
#' @param sortby rank the GESS result by one of the following scores: `WTCS`, `NCS`, `Tau`, `NCSct` or `NA` 
#' @param chunk_size size of chunk per processing
#' @return data.frame of GESS result, a list of drugs in reference database ranked by their similarity to query signature
#' @importFrom utils download.file
#' @importFrom R.utils gunzip
#' @import HDF5Array
#' @import SummarizedExperiment
#' @export
#' 
gess_lincs <- function(qSig, ES_NULL="Default", taurefList="Default", sortby="NCS", chunk_size=5000){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  #stopifnot(validObject(qSig))
  if(qSig@gess_method != "LINCS"){
    stop("The 'gess_method' slot of 'qSig' should be 'LINCS' if using 'gess_lincs' function")
  }
  upset <- qSig@qsig[[1]]
  downset <- qSig@qsig[[2]]
  se <- qSig@refdb 
  res <- lincsEnrich(se=se, upset=upset, downset=downset, ES_NULL=ES_NULL, taurefList=taurefList, sortby=sortby, chunk_size=chunk_size)
  x <- gessResult(result = as_tibble(res),
                  qsig = qSig@qsig,
                  gess_method = qSig@gess_method,
                  refdb = qSig@refdb)
  return(x)
}

randQuerySets <- function(id_names, N_queries, set_length=150) {
  randset_names <- paste0("randset_", sprintf("%09d", seq_len(N_queries)))
  rand_query_list <- sapply(randset_names, function(x) {
    id_list <- sample(id_names, 2 * set_length)
    split(id_list, rep(c("up", "down"), each=set_length))
  },
  simplify=FALSE
  )
}

#' ES_null Distribution with Random Queryies for Computing Nominal P-values for ES
#' @param se SummarizedExperiment object loaded from HDF5 backend LICNS GSE92742 database via `loadHDF5SummarizedExperiment(dbpath)` function
#' `dbpath` is the direcotry path to the `lincs42` database
#' @param N_queries number of random queries
#' @param dest_ES_NULL_path file path to the generated "ES_NULL.txt" file
#' @param partition partition node used to submit jobs on cluster
#' @import batchtools
#' @importFrom utils write.table
#' @return ES_NULL.tsv file
#' @export

randQueryES_slurm <- function(se, N_queries=1000, dest_ES_NULL_path, partition) {
  ## Create list of random queries
  idnames <- rownames(se)
  query_list <- randQuerySets(id_names=idnames, N_queries=N_queries, set_length=150)
  # saveRDS(query_list, query_list_path)
  ## Define submission function
  f <- function(x, query_list, se, dest_ES_NULL_path) {
    chunkno <- x 
    sz <- 10 # small enough to run within 10 hours 
    qlist <- split(query_list, ceiling(seq_along(names(query_list))/sz)) 
    qlistone <- qlist[[chunkno]] 
    for(i in seq_along(qlistone)) {
      esout <- lincsEnrich(se, upset=qlistone[[i]]$up, downset=qlistone[[i]]$down, sortby=NA, output="esonly", type=1)
      names(esout) <- colData(se)$pert_cell_factor
      wtcs = esout
      if(i==1) myMA <- matrix(, length(esout), length(qlistone), dimnames=list(names(wtcs), seq_along(qlistone)))
      myMA[,i] <- wtcs[rownames(myMA)]
      colnames(myMA)[i] <- names(qlistone[i])
    }
    return(myMA)
  }
  ## Split query into chunks, each chunk will be processed on cluster as one process
  sz <- 10 # small enough to use short queue 
  qlist <- split(query_list, ceiling(seq_along(names(query_list))/sz)) 
  dir = dirname(dest_ES_NULL_path)
  setwd(dir)
  if(! file.exists("slurm.tmpl")) download.file("https://goo.gl/tLMddb", "slurm.tmpl", quiet=TRUE)
  if(! file.exists(".batchtools.conf.R")) download.file("https://goo.gl/5HrYkE", ".batchtools.conf.R", quiet=TRUE)
  if(dir.exists("es_queries_reg")) unlink("es_queries_reg", recursive=TRUE)
  reg = makeRegistry(file.dir="es_queries_reg", conf.file=".batchtools.conf.R", packages="signatureSearch")
  ids = batchMap(f, x = seq(along=qlist), more.args = list(query_list=query_list, se=se, dest_ES_NULL_path=dest_ES_NULL_path))
  #testJob(id = 1)
  done <- submitJobs(ids, resources=list(walltime=43200, ncpus=1, memory=2048, partition=partition), reg=reg)
  print(waitForJobs())
  #getJobTable()
  print(getStatus())
  print(getErrorMessages())
  #showLog(1)
  #outlist <- batchtools::btlapply(seq(along=qlist), f, more.args = list(query_list=query_list, se=se, dest_ES_NULL_path=dest_ES_NULL_path), resources=list(walltime=600, ntasks=1, ncpus=1, memory=2048, partition=partition))
  # removeRegistry(wait=0, reg=reg)
  
  ## Collect results in frequency table with 3 diget accuracy
  if(waitForJobs()){
    esMA <- data.frame(WTCS=as.character(round(rev(seq(-1, 1, by=0.001)), 3)), Freq=0, stringsAsFactors=FALSE)
    for(i in ids$job.id) {
      mat <- loadResult(i)
      freq <- table(round(as.numeric(mat,3),3)) # processes entire data.frame
      freq <- freq[as.character(esMA[,1])]
      freq[is.na(freq)] <- 0
      esMA[,"Freq"] <- as.numeric(esMA[,"Freq"]) + as.numeric(freq)
      print(paste("Processed", i))
    }
    write.table(esMA, file=dest_ES_NULL_path, quote=FALSE, row.names=FALSE, sep="\t")
  } else {
    print("Not all the submitted jobs are successfully completed, please check!")
  }
}

#' ES_null Distribution with Random Queryies for Computing Nominal P-values for ES
#' @param se SummarizedExperiment object loaded from HDF5 backend LICNS GSE92742 database via `loadHDF5SummarizedExperiment(dbpath)` function
#' `dbpath` is the direcotry path to the `lincs42` database
#' @param N_queries number of random queries
#' @param dest_ES_NULL_path file path to the generated "ES_NULL.txt" file
#' @importFrom utils write.table
#' @return ES_NULL.txt file
#' @export

randQueryES_local <- function(se, N_queries=1000, dest_ES_NULL_path) {
  ## Create list of random queries
  idnames <- rownames(se)
  query_list <- randQuerySets(id_names=idnames, N_queries=N_queries, set_length=150)
  ## Define sapply function
  f <- function(x, query_list, se) {
      esout <- lincsEnrich(se, upset=query_list[[x]]$up, downset=query_list[[x]]$down, sortby=NA, output="esonly", type=1)
      names(esout) <- colData(se)$pert_cell_factor
      message("Random query ", sprintf("%04d", x), " has been searched against reference database")
      wtcs = esout
  }
  myMA <- sapply(seq(along=query_list), f, query_list, se)
  colnames(myMA) <- names(query_list)
  # write.table(as.data.frame(myMA), file=paste0(dirname(dest_ES_NULL_path), "/myMA.tsv"), col.names = NA, quote=FALSE, sep="\t")
  
  ## Collect results in frequency table with 3 diget accuracy
  esMA <- data.frame(WTCS=as.character(round(rev(seq(-1, 1, by=0.001)), 3)), Freq=0, stringsAsFactors=FALSE) 
  freq <- table(round(as.numeric(as.matrix(myMA),3),3)) # processes entire myMA data.frame
  freq <- freq[as.character(esMA[,1])]
  freq[is.na(freq)] <- 0
  esMA[,"Freq"] <- as.numeric(esMA[,"Freq"]) + as.numeric(freq)
  write.table(esMA, file=dest_ES_NULL_path, quote=FALSE, row.names=FALSE, sep="\t")
}

#' Create Query Reference DB for Tau Score Computation of lincsEnrich
#' @param se `SummarizedExperiment` object represents refernce database
#' @param Nup number of up-regulated genes subsetted for each signature in reference database `se` as query signature
#' @param Ndown number of down-regulated genes subsetted for each signature in reference database `se` as query signature
#' @param ES_NULL path to the ES_NULL file. ES_null distribution is generated with random queryies for computing nominal P-values for ES by 
#' `randQueryES_slurm` or `randQueryES_local` function. If `ES_NULL` is set as `Default`, it uses the ES_null distribution that we generated.
#' @param dest_path path to the destination file, including the file name "taurefList.rds"
#' @param partition partition node used to submit jobs on cluster
#' @return creat "taurefList.rds" file
#' @export
queryReferenceDB <- function(se, Nup=150, Ndown=150, ES_NULL="Default", dest_path, partition) {
  dir = dirname(dest_path)
  score_mat = assay(se)
  ## Create query list for all signatures in se
  query_list <- lapply(colnames(score_mat), function(x) {
    sigvec = sort(as.matrix(score_mat[,x]), decreasing = TRUE)
    list(upset=utils::head(names(sigvec), Nup), downset=utils::tail(names(sigvec), Ndown))
  })
  names(query_list) = colData(se)$pert_cell_factor
  ## Define submission function
  f <- function(x, se, query_list, ES_NULL, dest_path) {
    chunkno <- x 
    sz <- 10 # small enough to use short queue 
    qlist <- split(query_list, ceiling(seq_along(names(query_list))/sz))
    myMA <- matrix(, length(query_list), sz, dimnames=list(names(query_list), 1:sz))
    qlistone <- qlist[[chunkno]] 
    for(i in seq_along(qlistone)) {
      resultDF <- lincsEnrich(se, upset=qlistone[[i]]$up, downset=qlistone[[i]]$down, sortby=NA, output="no_tau", ES_NULL=ES_NULL, taurefList="Default")
      ncs <- resultDF$NCS
      mynames <- paste(resultDF$Pert, resultDF$Type, sep="__")
      mynames <- gsub("__up|__down", "", mynames)
      names(ncs) <- mynames
      myMA[,i] <- ncs[rownames(myMA)]
      colnames(myMA)[i] <- names(qlistone[i])
    }
    myMA <- myMA[, seq_along(qlistone)] # Only relevant for last entry that may not have as many colums as sz
    return(myMA)
  }
  ## Split query into chunks, each chunk will be processed on cluster as one process
  sz <- 10 # small enough to use short queue 
  qlist <- split(query_list, ceiling(seq_along(names(query_list))/sz)) 
  # qlist <-  qlist[1:200] # test
  # qlist <-  qlist[201:length(qlist)] # test
  dir = dirname(dest_path)
  setwd(dir)
  if(! file.exists("slurm.tmpl")) download.file("https://goo.gl/tLMddb", "slurm.tmpl", quiet=TRUE)
  if(! file.exists(".batchtools.conf.R")) download.file("https://goo.gl/5HrYkE", ".batchtools.conf.R", quiet=TRUE)
  if(dir.exists("tau_queries_reg")) unlink("tau_queries_reg", recursive=TRUE)
  reg = makeRegistry(file.dir="tau_queries_reg", conf.file=".batchtools.conf.R", packages="signatureSearch")
  ids = batchMap(f, x = seq(along=qlist), more.args = list(se=se, query_list=query_list, ES_NULL=ES_NULL, dest_path=dest_path))
  #testJob(id = 1)
  done <- submitJobs(ids, resources=list(walltime=43200, ncpus=1, memory=2048, partition=partition), reg=reg)
  print(waitForJobs())
  #getJobTable()
  print(getStatus())
  print(getErrorMessages())
  #showLog(1)
  # removeRegistry(wait=0, reg=reg)
  
  ## Inspect result and resubmit jobs for missing and empty files
  #dir = dirname(dest_path)
  #fileDF <- file.info(list.files(paste0(dir, "/tau_queries"), pattern="result_*", full.names=TRUE))
  #index_empty <- as.numeric(gsub(".*_", "", row.names(fileDF[fileDF$size==0, ])))
  #qlist <- split(query_list, ceiling(seq_along(names(query_list))/sz)) 
  #index_all_files <- seq_along(qlist)
  #index_exist <- as.numeric(gsub(".*_", "", row.names(fileDF)))
  #index_missing <- index_all_files[!index_all_files %in% index_exist]
  #index_repeat <- unique(sort(c(index_empty, index_missing)))
  #if(length(index_repeat)!=0) outlist <- bplapply(index_repeat, f)
  
  ## Organize results in list where each component contains data.frame with query results from one cell type
  if(waitForJobs){
    pathDF <- data.frame(query=names(query_list), path=rep(seq_along(qlist), each=sz))
    pathDF <- data.frame(pathDF, target=gsub("^.*?__", "", pathDF$query))
    pathList <- split(as.character(pathDF$path), factor(pathDF$target))
    pathList <- sapply(pathList, unique) # eliminates unnecessary/duplicated imports of files
    taurefList <- rep(NA, length(pathList)); names(taurefList) <- names(pathList)
    taurefList <- as.list(taurefList) 
    for(celltype in names(pathList)) {
      for(i in seq_along(pathList[[celltype]])) {
        tmpDF <- loadResult(as.numeric(pathList[[celltype]][i]))
        tmpDF <- round(tmpDF, 2) # Reduces data size
        colindex <- gsub("^.*?__", "", colnames(tmpDF)) %in% celltype
        tmpDF <- tmpDF[, colindex, drop=FALSE] 
        if(i==1) { 
          rowindex <- gsub("^.*?__", "", rownames(tmpDF)) %in% celltype
          containerDF <- tmpDF[rowindex, , drop=FALSE]
        } else {
          containerDF <- cbind(containerDF, tmpDF[rownames(containerDF),])
        }
        print(paste("Finished", i, "of", length(pathList[[celltype]]), "results collected in", ncol(containerDF), "columns."))
      }
      taurefList[[celltype]] <- containerDF
      rm(containerDF); gc()
    }
    saveRDS(taurefList, file=dest_path)
  } else {
    print("Not all the submitted jobs are successfully completed, please check!")
  }
}
