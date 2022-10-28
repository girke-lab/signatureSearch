#' @title GESS Methods
#' @rdname gess
#' @description 
#' LINCS search method implements the Gene Expression Signature Search (GESS) from 
#' Subramanian et al, 2017, here referred to as LINCS. The method uses as 
#' query the two label sets of the most up- and down-regulated genes from a 
#' genome-wide expression experiment, while the reference database is composed 
#' of differential gene expression values (e.g. LFC or z-scores). Note, the 
#' related CMAP method uses here ranks instead.
#' @details 
#' Subramanian et al. (2017) introduced a more complex GESS algorithm, 
#' here referred to as LINCS. While related to CMAP, there are several important
#' differences among the two approaches. First, LINCS weights the query genes 
#' based on the corresponding differential expression scores of the GESs in the 
#' reference database (e.g. LFC or z-scores). Thus, the reference database used 
#' by LINCS needs to store the actual score values rather than their ranks.
#' Another relevant difference is that the LINCS algorithm uses a bi-directional
#' weighted Kolmogorov-Smirnov enrichment statistic (ES) as similarity metric.
#' @section Column description:
#' Descriptions of the columns in GESS result tables.
#' \itemize{
#'   \item pert: character, perturbagen (e.g. drugs) in the reference 
#'   database. The treatment/column names of the reference database are 
#'   organized as \code{pert__cell__trt_cp} format. The \code{pert} column in
#'   GESS result table contains what stored under the \code{pert} slot of the 
#'   column names. 
#'   \item cell: character, acronym of cell type
#'   \item type: character, perturbation type. In the CMAP/LINCS 
#'   databases provided by \code{signatureSearchData}, the perturbation types
#'   are currently treatments with drug-like compounds (trt_cp). If required,
#'   users can build custom signature database with other types of
#'   perturbagens (e.g., gene knockdown or over-expression events) with the 
#'   provided \code{\link{build_custom_db}} function.
#'   \item trend: character, up or down when the reference signature is 
#'   positively or negatively connected with the query signature, 
#'   respectively.
#'   \item N_upset: integer, number of genes in the query up set
#'   \item N_downset: integer, number of genes in the query down set
#'   \item WTCS: Weighted Connectivity Score, a bi-directional Enrichment 
#'   Score for an up/down query set. If the ES values of an up set and a down 
#'   set are of different signs, then WTCS is (ESup-ESdown)/2, otherwise, 
#'   it is 0. WTCS values range from -1 to 1. They are positive or negative 
#'   for signatures that are positively or inversely related, respectively, 
#'   and close to zero for signatures that are unrelated.
#'   \item WTCS_Pval: Nominal p-value of WTCS computed by comparing WTCS 
#'   against a null distribution of WTCS values obtained from a large number
#'   of random queries (e.g. 1000).
#'   \item WTCS_FDR: False discovery rate of WTCS_Pval.
#'   \item NCS: Normalized Connectivity Score. To make connectivity scores 
#'   comparable across cell types and perturbation types, 
#'   the scores are normalized. Given a vector of WTCS 
#'   values w resulting from a query, the values are normalized within each 
#'   cell line c and perturbagen type t to obtain NCS by dividing the WTCS 
#'   value with the signed mean of the WTCS values within 
#'   the subset of the signatures in the reference database corresponding to c 
#'   and t.
#'   \item Tau: Enrichment score standardized for a given database. 
#'   The Tau score compares an observed NCS to a large set of NCS values that 
#'   have been pre-computed for a specific reference database. The query results
#'   are scored with Tau as a standardized measure ranging from 100 to -100. 
#'   A Tau of 90 indicates that only 10% of reference perturbations exhibit 
#'   stronger connectivity to the query. This way one can make more meaningful 
#'   comparisons across query results. 
#'   
#'   Note, there are NAs in the Tau score column, the reason is that the number 
#'   of signatures in \emph{Qref} that match the cell line of signature \emph{r}
#'   (the \code{TauRefSize} column in the GESS result) is less than 500, 
#'   Tau will be set as NA since it is redeemed as there are not large enough 
#'   samples for computing meaningful Tau scores.
#'   
#'   \item TauRefSize: Size of reference perturbations for computing Tau.
#'   \item NCSct: NCS summarized across cell types. Given a vector of NCS values
#'   for perturbagen p, relative to query q, across all cell lines c in which p 
#'   was profiled, a cell-summarized connectivity score is obtained using a 
#'   maximum quantile statistic. It compares the 67 and 33 quantiles of 
#'   NCSp,c and retains whichever is of higher absolute magnitude.
#'   \item cor_score: Correlation coefficient based on the method defined in 
#'   the \code{gess_cor} function.
#'   \item raw_score: bi-directional enrichment score (Kolmogorov-Smirnov 
#'   statistic) of up and down set in the query signature
#'   \item scaled_score: raw_score scaled to values from 1 to -1 by 
#'   dividing the positive and negative scores with the maximum positive score 
#'   and the absolute value of the minimum negative score, respectively.
#'   \item effect: Scaled bi-directional enrichment score corresponding to 
#'   the scaled_score under the CMAP result.
#'   \item nSet: number of genes in the GES in the reference
#'   database (gene sets) after setting the higher and lower cutoff.
#'   \item nFound: number of genes in the GESs of the reference
#'   database (gene sets) that are also present in the query GES.
#'   \item signed: whether gene sets in the reference database have signs, 
#'   representing up and down regulated genes when computing scores. 
#'   \item pval: p-value of the Fisher's exact test.
#'   \item padj: p-value adjusted for multiple hypothesis testing using
#'   R's p.adjust function with the Benjamini & Hochberg (BH) method. 
#'   \item effect: z-score based on the standard normal distribution.
#'   \item LOR: Log Odds Ratio.
#'   \item t_gn_sym: character, symbol of the gene encoding the
#'   corresponding drug target protein
#'   \item MOAss: character, compound MOA annotation from \code{signatureSearch}
#'   package
#'   \item PCIDss: character, compound PubChem CID annotation from 
#'   \code{signatureSearch} package
#' }
#' 
#' @param qSig \code{\link{qSig}} object defining the query signature including
#' the GESS method (should be 'LINCS') and the path to the reference database.
#' For details see help of \code{qSig} and \code{qSig-class}.
#' @param tau TRUE or FALSE, whether to compute the tau score. Note, TRUE is 
#' only meaningful when the full LINCS database is searched, since accurate Tau 
#' score calculation depends on the usage of the exact same database their 
#' background values are based on.
#' @param sortby sort the GESS result table based on one of the following 
#' statistics: `WTCS`, `NCS`, `Tau`, `NCSct` or `NA` 
#' @param chunk_size number of database entries to process per iteration to 
#' limit memory usage of search.
#' @param ref_trts character vector. If users want to search against a subset 
#' of the reference database, they could set ref_trts as a character vector 
#' representing column names (treatments) of the subsetted refdb. 
#' @param workers integer(1) number of workers for searching the reference
#' database parallelly, default is 1.
#' @param method One of 'spearman' (default), 'kendall', or 'pearson',
#' indicating which correlation coefficient to use.
#' 
#' @param higher The 'upper' threshold. If not 'NULL', genes with a score 
#' larger than or equal to 'higher' will be included in the gene set with 
#' sign +1. At least one of 'lower' and 'higher' must be specified. 
#' 
#' \code{higher} argument need to be set as \code{1} if the \code{refdb} in 
#' \code{qSig} is path to the HDF5 file that were converted from the gmt file.
#' 
#' @param lower The lower threshold. If not 'NULL', genes with a score smaller 
#' than or equal 'lower' will be included in the gene set with sign -1. 
#' At least one of 'lower' and 'higher' must be specified.
#' 
#' \code{lower} argument need to be set as \code{NULL} if the \code{refdb} in 
#' \code{qSig} is path to the HDF5 file that were converted from the gmt file.
#' 
#' @param padj numeric(1), cutoff of adjusted p-value or false discovery rate 
#' (FDR) of defining DEGs that is less than or equal to 'padj'. The 'padj' 
#' argument is valid only if the reference HDF5 file contains the p-value 
#' matrix stored in the dataset named as 'padj'. 
#' @param addAnnotations Logical value.  If \code{TRUE} adds drug annotations to results. 
#' @inheritParams addGESSannot
#' 
#' @return \code{\link{gessResult}} object, the result table contains the 
#' search results for each perturbagen in the reference database ranked by 
#' their signature similarity to the query.
#' @import SummarizedExperiment
#' @seealso \code{\link{qSig}}, \code{\link{gessResult}}, 
#'          \code{\link{addGESSannot}}
#' @references For detailed description of the LINCS method and scores, 
#' please refer to: Subramanian, A., Narayan, R., Corsello, S. M., Peck, D. D., 
#' Natoli, T. E., Lu, X., Golub, T. R. (2017). A Next Generation 
#' Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell, 
#' 171 (6), 1437-1452.e17. URL: https://doi.org/10.1016/j.cell.2017.10.049
#' 
#' For detailed description of the CMap method, please refer to: 
#' Lamb, J., Crawford, E. D., Peck, D., Modell, J. W., Blat, I. C., 
#' Wrobel, M. J., Golub, T. R. (2006). The Connectivity Map: 
#' using gene-expression signatures to connect small molecules, genes, and 
#' disease. Science, 313 (5795), 1929-1935. 
#' URL: https://doi.org/10.1126/science.1132939
#' 
#' Sandmann, T., Kummerfeld, S. K., Gentleman, R., & Bourgon, R. 
#' (2014). gCMAP: user-friendly connectivity mapping with R. Bioinformatics , 
#' 30 (1), 127-128. URL: https://doi.org/10.1093/bioinformatics/btt592
#' 
#' Graham J. G. Upton. 1992. Fisher's Exact Test. J. R. Stat. Soc. Ser. A 
#' Stat. Soc. 155 (3). [Wiley, Royal Statistical Society]: 395-402. 
#' URL: http://www.jstor.org/stable/2982890
#' @examples 
#' 
#' ############### LINCS method #############
#' # qsig_lincs <- qSig(query=list(
#' #                      upset=c("230", "5357", "2015", "2542", "1759"), 
#' #                      downset=c("22864", "9338", "54793", "10384", "27000")), 
#' #                    gess_method="LINCS", refdb=db_path)
#' # lincs <- gess_lincs(qsig_lincs, sortby="NCS", tau=FALSE)
#' # result(lincs)
#' @export
gess_lincs <- function(qSig, tau=FALSE, sortby="NCS", 
                       chunk_size=5000, ref_trts=NULL, workers=1,
                       cmp_annot_tb=NULL, by="pert", cmp_name_col="pert", addAnnotations = TRUE){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  if(gm(qSig) != "LINCS"){
    stop(paste("The 'gess_method' slot of 'qSig' should be 'LINCS'",
               "if using 'gess_lincs' function"))
  }
  upset <- qr(qSig)$upset
  downset <- qr(qSig)$downset
  db_path <- determine_refdb(refdb(qSig))
  res <- lincsEnrich(db_path, upset=upset, downset=downset, 
                     tau=tau, sortby=sortby, chunk_size=chunk_size, 
                     ref_trts=ref_trts, workers=workers)
  # add compound annotations
  if(addAnnotations == TRUE){
  res <- addGESSannot(res, refdb(qSig), cmp_annot_tb =cmp_annot_tb[,!colnames(cmp_annot_tb) %in% "t_gn_sym"], by, cmp_name_col)
  } else {
    res <- as_tibble(res)  
  }
  
  x <- gessResult(result = res,
                  query = qr(qSig),
                  gess_method = gm(qSig),
                  refdb = refdb(qSig))
  return(x)
}


lincsEnrich <- function(db_path, upset, downset, sortby="NCS", type=1, 
                        output="all", tau=FALSE, minTauRefSize=500, 
                        chunk_size=5000, ref_trts=NULL, workers=4) {
    mycolnames <- c("WTCS", "NCS", "Tau", "NCSct", "N_upset", "N_downset", NA)
    if(!any(mycolnames %in% sortby)) 
        stop("Unsupported value assinged to sortby.")
    
    ## calculate ESout of query to blocks (e.g., 5000 columns) of full refdb
    full_mat <- HDF5Array(db_path, "assay")
    rownames(full_mat) <- as.character(HDF5Array(db_path, "rownames"))
    colnames(full_mat) <- as.character(HDF5Array(db_path, "colnames"))
    
    if(! is.null(ref_trts)){
        trts_valid <- trts_check(ref_trts, colnames(full_mat))
        full_mat <- full_mat[, trts_valid]
    }
    
    full_dim <- dim(full_mat)
    full_grid <- colAutoGrid(full_mat, ncol=min(chunk_size, ncol(full_mat)))
    ### The blocks in 'full_grid' are made of full columns 
    nblock <- length(full_grid) 
    
    ESout <- unlist(bplapply(seq_len(nblock), function(b){
      ref_block <- read_block(full_mat, full_grid[[as.integer(b)]])
      mat <- ref_block
      ## Run .enrichScore on upset and downset
      ## When both upset and downset are provided 
      if(length(upset)>0 & length(downset)>0) {
        ESup <- apply(mat, 2, function(x) 
            .enrichScore(sigvec=sort(x, decreasing = TRUE), 
                         Q=upset, type=type))
        ESdown <- apply(mat, 2, function(x) 
            .enrichScore(sigvec=sort(x, decreasing = TRUE),
                         Q=downset, type=type))
        ESout1 <- ifelse(sign(ESup) != sign(ESdown), (ESup - ESdown)/2, 0)
        ## When only upset is provided
      } else if(length(upset)>0 & length(downset)==0) {
        ESup <- apply(mat, 2, function(x) 
          .enrichScore(sigvec=sort(x, decreasing = TRUE), 
                       Q=upset, type=type))
        ESout1 <- ESup
        ## When only downset is provided
      } else if(length(upset)==0 & length(downset)>0) {
        ESdown <- apply(mat, 2, function(x) 
          .enrichScore(sigvec=sort(x, decreasing = TRUE), 
                       Q=downset, type=type))
        ESout1 <- -ESdown
        ## When none are provided (excluded by input validity check already)
      }
      }, BPPARAM=MulticoreParam(workers=workers)))
    
    # ## Read in matrix in h5 file by chunks
    # mat_dim <- getH5dim(db_path)
    # mat_nrow <- mat_dim[1]
    # mat_ncol <- mat_dim[2]
    # ceil <- ceiling(mat_ncol/chunk_size)
    # ESout <- NULL
    # for(i in seq_len(ceil)){
    #     mat <- readHDF5mat(db_path,
    #               colindex=(chunk_size*(i-1)+1):min(chunk_size*i, mat_ncol))
    #     ## Run .enrichScore on upset and downset
    #     ## When both upset and downset are provided 
    #     if(length(upset)>0 & length(downset)>0) {
    #         ESup <- apply(mat, 2, function(x) 
    #             .enrichScore(sigvec=sort(x, decreasing = TRUE), 
    #                          Q=upset, type=type))
    #         ESdown <- apply(mat, 2, function(x) 
    #             .enrichScore(sigvec=sort(x, decreasing = TRUE),
    #                          Q=downset, type=type)) 
    #         ESout1 <- ifelse(sign(ESup) != sign(ESdown), (ESup - ESdown)/2, 0)
    #         ## When only upset is provided
    #     } else if(length(upset)>0 & length(downset)==0) {
    #         ESup <- apply(mat, 2, function(x) 
    #             .enrichScore(sigvec=sort(x, decreasing = TRUE), 
    #                          Q=upset, type=type))
    #         ESout1 <- ESup
    #         ## When only downset is provided
    #     } else if(length(upset)==0 & length(downset)>0) {
    #         ESdown <- apply(mat, 2, function(x) 
    #             .enrichScore(sigvec=sort(x, decreasing = TRUE), 
    #                          Q=downset, type=type))
    #         ESout1 <- -ESdown
    #         ## When none are provided (excluded by input validity check already)
    #     }
    #     ESout <- c(ESout, ESout1)
    # }
    
    ## Assmble output 
    if(output=="esonly") {
        return(ESout)
    }
    if(output=="all") {
        resultDF <- .lincsScores(esout=ESout, upset=upset, downset=downset, 
                                 minTauRefSize=minTauRefSize, tau=tau)
    }
    if(!is.na(sortby)) {
        resultDF <- resultDF[order(abs(resultDF[,sortby]), decreasing=TRUE), ]
    } else {
        resultDF <- resultDF
    }
    row.names(resultDF) <- NULL
    return(resultDF)
}

#' @importFrom utils read.delim
#' @importFrom stats quantile
.lincsScores <- function(esout, upset, downset, minTauRefSize, tau=FALSE) {
    ## P-value and FDR for WTCS based on ESnull from random queries where 
    ## p-value = sum(ESrand > ES_obs)/Nrand
    
    # download ES_NULL.txt from AnnotationHub
    WTCSnull <- validLoad("EH3234")
    WTCSnull[WTCSnull[, "Freq"]==0, "Freq"] <- 1 
    # Add pseudo count of 1 where Freq is zero 
    myrounding <- max(nchar(as.character(WTCSnull[,"WTCS"]))) - 3 
    # Three because of dot and minus sign
    es_round <- round(as.numeric(esout), myrounding) 
    # Assures same rounding used for WTCSnull computation
    WTCS_pval <- vapply(es_round, function(x) {
      sum(WTCSnull[abs(WTCSnull[,"WTCS"]) > abs(x),"Freq"])/sum(WTCSnull[,"Freq"])
      }, FUN.VALUE = numeric(1))
    WTCS_fdr <- p.adjust(WTCS_pval, "fdr")
    ## Normalized connectivity score (NCS)
    grouping <- paste(gsub("^.*?__", "", names(esout)), 
                      as.character(ifelse(esout > 0, "up", "down")), sep="__")
    es_na <- as.numeric(esout)
    es_na[es_na == 0] <- NA 
    # eliminates zeros from mean calculation; zeros have high impact on NCS 
    # values due to their high frequency
    groupmean <- tapply(es_na, grouping, mean, na.rm=TRUE)
    groupmean[is.na(groupmean)] <- 0 
    # In case groups contain only zeros, NA/NaN values are introduced in mean 
    # calculation which are reset to zeros here
    groupmean[groupmean==0] <- 10^-12 
    # Set zeros (can be from non NAs) to small value to avoid devision by zero 
    ncs <- as.numeric(esout) / abs(groupmean[grouping]) 
    # without abs() sign of neg values would switch to pos
    ## Tau calculation requires reference NCS lookup DB
    ## performs: sign(ncs_query) * 100/N sum(abs(ncs_ref) < abs(ncs_query))
    if(tau){
      # download taurefList.rds
      taurefList9264 <- validLoad("EH3233")
      
      ncs_query <- ncs; names(ncs_query) <- names(esout)
      queryDB_refDB_match <- 
          unique(unlist(lapply(taurefList9264, rownames))) %in% names(ncs_query)
      if(!all(queryDB_refDB_match)) warning(
          paste0("QueryDB and tauRefDB differ by ", 
          round(100 * sum(!queryDB_refDB_match)/length(queryDB_refDB_match),1), 
          "% of their entries.",
          " Accurate tau computation requires close to 0% divergence. \n"))
      ncs_query_list <- split(ncs_query, 
                              factor(gsub("^.*?__", "", names(ncs_query))))
      tau_score <- lapply(names(ncs_query_list), function(x) {
        tmpDF <- taurefList9264[[x]]
        ncs_query_match <- names(ncs_query_list[[x]])[names(ncs_query_list[[x]]) 
                                                      %in% rownames(tmpDF)]
        
        if(length(ncs_query_match)>0) {
          tmpDF <- tmpDF[ncs_query_match, , drop=FALSE]
          #### subset to the same length ##### rounded as in ref db
          sign(ncs_query_list[[x]][names(ncs_query_list[[x]]) %in% rownames(tmpDF)]) * 100/ncol(tmpDF) * 
            rowSums(abs(tmpDF) < abs(round(ncs_query_list[[x]][names(ncs_query_list[[x]]) %in% rownames(tmpDF)], 2))) 
        } else {
          NULL
        }
      })
      tau_score <- unlist(tau_score)
      tau_score <- tau_score[names(ncs_query)]
      tauRefSize <- vapply(taurefList9264, ncol, 
                  FUN.VALUE = integer(1))[gsub("^.*?__", "", names(tau_score))]
      tau_score[tauRefSize < minTauRefSize] <- NA
      ## Add by YD 
      rm(taurefList9264); gc()
    }
    ## Summary across cell lines (NCSct)
    ctgrouping <- gsub("__.*__", "__", names(esout))
    qmax <- tapply(ncs, ctgrouping, function(x) { 
      q <- quantile(x, probs=c(0.33, 0.67))
      ifelse(abs(q[2]) >= abs(q[1]), q[2], q[1])
    })
    qmax <- qmax[ctgrouping]
    ## Organize result in data.frame
    new <- as.data.frame(t(vapply(seq_along(esout), function(i)
      unlist(strsplit(as.character(names(esout)[i]), "__")),
      FUN.VALUE = character(3))), stringsAsFactors=FALSE)
    colnames(new) <- c("pert", "cell", "type")
    
    if(tau){
      resultDF <- data.frame(
        new, 
        trend = as.character(ifelse(esout > 0, "up", "down")),
        WTCS = as.numeric(esout), 
        WTCS_Pval = WTCS_pval,
        WTCS_FDR = WTCS_fdr, 
        NCS = ncs, 
        Tau = tau_score,
        TauRefSize=tauRefSize, 
        NCSct = qmax, 
        N_upset = length(upset),
        N_downset = length(downset), stringsAsFactors = FALSE)
    } else {
      resultDF <- data.frame(
        new, 
        trend = as.character(ifelse(esout > 0, "up", "down")),
        WTCS = as.numeric(esout), 
        WTCS_Pval = WTCS_pval,
        WTCS_FDR = WTCS_fdr, 
        NCS = ncs, 
        NCSct = qmax, 
        N_upset = length(upset),
        N_downset = length(downset), stringsAsFactors = FALSE)
    }
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
    if(NR == 0) return(0)
    ESvec <- cumsum((hit_index * R * 1/NR) - (miss_index * 1/Ns)) 
    ES <- ESvec[which.max(abs(ESvec))]
    return(ES)
}

#' Function computes null distribution of Weighted Connectivity Scores (WTCS) 
#' used by the LINCS GESS method for computing nominal P-values.
#' 
#' @title Generate WTCS Null Distribution with Random Queries
#' @param h5file character(1), path to the HDF5 file representing the
#' reference database
#' @param N_queries number of random queries
#' @param dest path to the output file (e.g. "ES_NULL.txt")
#' @return File with path assigned to \code{dest}
#' @importFrom utils write.table
#' @examples 
#' db_path = system.file("extdata", "sample_db.h5", package="signatureSearch")
#' rand_query_ES(h5file=db_path, N_queries=5, dest="ES_NULL.txt")
#' unlink("ES_NULL.txt")
#' @seealso \code{\link{gess_lincs}}
#' @references 
#' Subramanian, A., Narayan, R., Corsello, S. M., Peck, D. D., Natoli, T. E., 
#' Lu, X., Golub, T. R. (2017). A Next Generation Connectivity Map: L1000 
#' Platform and the First 1,000,000 Profiles. Cell, 171 (6), 1437-1452.e17. 
#' URL: https://doi.org/10.1016/j.cell.2017.10.049
#' @export

rand_query_ES <- function(h5file, N_queries=1000, dest) {
  ## Create list of random queries
  idnames <- drop(h5read(h5file, "rownames"))
  query_list <- randQuerySets(id_names=idnames, N_queries=N_queries, 
                              set_length=150)
  ## Define vapply function
  f <- function(x, query_list, h5file) {
      esout <- lincsEnrich(h5file, upset=query_list[[x]]$up, 
                           downset=query_list[[x]]$down, sortby=NA, 
                           output="esonly", type=1)
      names(esout) <- drop(h5read(h5file, "colnames"))
      # message("Random query ", sprintf("%04d", x), 
      #         " has been searched against reference database")
      wtcs <- esout
  }
  myMA <- vapply(seq(along=query_list), f, query_list, h5file, 
                 FUN.VALUE=double(length(drop(h5read(h5file, "colnames")))))
  colnames(myMA) <- names(query_list)
  ## Collect results in frequency table with 3 diget accuracy
  esMA <- data.frame(WTCS=as.character(round(rev(seq(-1, 1, by=0.001)), 3)), 
                     Freq=0, stringsAsFactors=FALSE) 
  freq <- table(round(as.numeric(as.matrix(myMA),3),3)) 
  ## processes entire myMA data.frame
  freq <- freq[as.character(esMA[,1])]
  freq[is.na(freq)] <- 0
  esMA[,"Freq"] <- as.numeric(esMA[,"Freq"]) + as.numeric(freq)
  write.table(esMA, file=dest, quote=FALSE, 
              row.names=FALSE, sep="\t")
}

randQuerySets <- function(id_names, N_queries, set_length=150) {
    randset_names <- paste0("randset_", sprintf("%09d", seq_len(N_queries)))
    rand_query_list <- lapply(randset_names, function(x) {
        id_list <- sample(id_names, 2 * set_length)
        split(id_list, rep(c("up", "down"), each=set_length))
    })
    names(rand_query_list) <- randset_names
    return(rand_query_list)
}

