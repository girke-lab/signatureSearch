#' @rdname gess
#' @description 
#' In its iterative form, Fisher's exact test (Upton, 1992) can be used as Gene 
#' Expression Signature (GES) Search to scan GES databases for entries that 
#' are similar to a query GES.
#' @details 
#' When using the Fisher's exact test (Upton, 1992) as GES Search (GESS) method,
#' both the query and the database are composed of gene label sets, such as 
#' DEG sets.
#' @importMethodsFrom GSEABase GeneSet
#' @importMethodsFrom GSEABase GeneSetCollection
#' @importMethodsFrom GSEABase incidence
#' @importFrom GSEABase EntrezIdentifier
#' @importFrom stats dhyper
#' @importFrom stats qnorm
#' @importFrom tibble tibble
#' @examples 
#' 
#' ############## Fisher's Exact Test ##########
#' # qsig_fisher <- qSig(query=query_mat, gess_method="Fisher", refdb=db_path)
#' # fisher <- gess_fisher(qSig=qsig_fisher, higher=1, lower=-1)
#' # result(fisher)
#' @export
#' 
gess_fisher <- function(qSig, higher=NULL, lower=NULL, padj=NULL,
                        chunk_size=5000, ref_trts=NULL, workers=1,
                        cmp_annot_tb=NULL, by="pert", cmp_name_col="pert",
                        addAnnotations = TRUE){
  if(!is(qSig, "qSig")) stop("The 'qSig' should be an object of 'qSig' class")
  #stopifnot(validObject(qSig))
  if(gm(qSig) != "Fisher"){
    stop("The 'gess_method' slot of 'qSig' should be 'Fisher' 
         if using 'gess_fisher' function")
  }
  if(is.list(qr(qSig))){
    query <- GeneSet(unique(unlist(qr(qSig))))
    query <- t(incidence(GeneSetCollection(query)))
  } else {
    query <- induceSMat(qr(qSig), higher=higher, lower=lower)
  }
  db_path <- determine_refdb(refdb(qSig))

  ## calculate fisher score of query to blocks (5000 columns) of full refdb
  full_mat <- HDF5Array(db_path, "assay")
  rownames(full_mat) <- as.character(HDF5Array(db_path, "rownames"))
  colnames(full_mat) <- as.character(HDF5Array(db_path, "colnames"))
  
  if(! is.null(ref_trts)){
    trts_valid <- trts_check(ref_trts, colnames(full_mat))
    full_mat <- full_mat[, trts_valid]
  }
  
  full_dim <- dim(full_mat)
  full_grid <- colAutoGrid(full_mat, ncol=min(chunk_size, ncol(full_mat)))
  
  if(! is.null(padj)){
      if(! 'padj' %in% h5ls(db_path)$name){
          stop("The 'refdb' need to be an hdf5 file that contains 'padj' dataset!")
      }
      full_pmat <- HDF5Array(db_path, name="padj")
      rownames(full_pmat) <- as.character(HDF5Array(db_path, name="rownames"))
      colnames(full_pmat) <- as.character(HDF5Array(db_path, name="colnames"))
  }
  
  ### The blocks in 'full_grid' are made of full columns 
  nblock <- length(full_grid) 
  resultDF <- bplapply(seq_len(nblock), function(b){
    ref_block <- read_block(full_mat, full_grid[[as.integer(b)]])
    if(! is.null(padj)){
        pmat <- read_block(full_pmat, full_grid[[as.integer(b)]])
    } else {
        pmat = NULL
    }
    cmap <- induceSMat(ref_block, higher=higher, lower=lower,
                       padj=padj, pmat=pmat)
    universe <- rownames(cmap)
    c <- fs(query=query, sets=cmap, universe = universe)
    return(data.frame(c))}, BPPARAM = MulticoreParam(workers = workers))
  resultDF <- do.call(rbind, resultDF)

  # else {
  #   gene_sets <- suppressWarnings(read_gmt(db_path))
  #   # remove invalid gene sets that have 0 length or names are NAs
  #   gene_sets <- gene_sets[sapply(gene_sets, length)>0 & !is.na(names(gene_sets))]
  #   # transform gene sets to sparseMatrix
  #   gsc <- GeneSetCollection(mapply(function(geneIds, setId) {
  #     GeneSet(geneIds, geneIdType=EntrezIdentifier(),
  #             setName=setId)
  #   }, gene_sets, names(gene_sets)))
  #   cmapData <- incidence(gsc)
  #   cmapData <- Matrix::t(cmapData)
  #   #c <- fs(query=query, sets=cmapData, universe=rownames(cmapData))
  #   
  #   # Search the refdb by using multiple cores/workers
  #   ceil <- ceiling(ncol(cmapData)/chunk_size)
  #   res_list <- bplapply(seq_len(ceil), function(i){
  #     dgcMat <- cmapData[,(chunk_size*(i-1)+1):min(chunk_size*i, ncol(cmapData))]
  #     c <- fs(query=query, sets=dgcMat, universe=rownames(dgcMat))
  #     return(data.frame(c))
  #   }, BPPARAM = MulticoreParam(workers=workers))
  #   resultDF <- do.call(rbind, res_list)
  # }
  
  # mat_dim <- getH5dim(db_path)
  # mat_ncol <- mat_dim[2]
  # ceil <- ceiling(mat_ncol/chunk_size)
  # fs <- function(i){
  #     mat <- readHDF5mat(db_path,
  #                colindex=(chunk_size*(i-1)+1):min(chunk_size*i, mat_ncol))
  #     cmap <- induceCMAPCollection(mat, higher=higher, lower=lower)
  #     universe <- featureNames(cmap)
  #     c <- fisher_score(query=query, sets=cmap, universe = universe)
  #     cmapTable(c)
  # }
  # resultDF <- do.call(rbind, lapply(seq_len(ceil), fs))
  
  resultDF <- resultDF[order(resultDF$padj), ]
  row.names(resultDF) <- NULL
  
  if(addAnnotations == TRUE){
  res <- sep_pcf(resultDF)
  # add compound annotations
  res <- addGESSannot(res, refdb(qSig), cmp_annot_tb = cmp_annot_tb[,!colnames(cmp_annot_tb) %in% "t_gn_sym"], by, cmp_name_col)
  } else {
    res <- tibble:::tibble(resultDF)
    colnames(res)[1] <- "pert"
  }
  x <- gessResult(result = res,
                  query = qr(qSig),
                  gess_method = gm(qSig),
                  refdb = refdb(qSig))
  return(x)
}

#' @importFrom Matrix crossprod
#' @importFrom Matrix colSums
fs <- function(query, sets, universe) {
    query <- query[intersect(rownames(query), universe), , drop=FALSE]
    sets <- sets[intersect(rownames(sets), universe), ]
    common.genes <- intersect(rownames(query), rownames(sets))
    
    if( length( common.genes) == 0 ) {
      stop( "None of the query gene identifiers could be found in the 
            reference dataset.",
            call. = FALSE)
    }
    
    query.common <- abs(matrix(query[common.genes,], nrow=length(common.genes)))
    sets.common <- abs(matrix(sets[common.genes,], nrow=length(common.genes)))
    coincidence <- Matrix::crossprod(query.common, sets.common) 
    
    ## 2x2 table
    ##              query  1         query  0
    ##   sets 1  query.and.sets  sets.not.query    || sets.all
    ##   sets 0  query.not.sets      neither 
    ##   ========================================
    ##            query.all        query.colsum
    
    query.and.sets <- t(matrix(coincidence, ncol=ncol(sets), 
                        dimnames=list(colnames(query), colnames(sets))))
    
    query.all <- matrix(rep(colSums(abs(query)), ncol(sets)),
                         nrow=ncol(sets), ncol=ncol(query), byrow=TRUE, 
                         dimnames=dimnames(query.and.sets))
    
    sets.all <- matrix(rep(colSums(abs(sets)), ncol(query)),
                       nrow=ncol(sets), ncol=ncol(query), byrow=FALSE, 
                       dimnames=dimnames(query.and.sets))
    
    neither <- length(universe) - sets.all - query.all + query.and.sets
    
    sets.not.query <- length(universe) - query.all - neither
    
    query.not.sets <- query.all - query.and.sets
    
    query.colsum <- sets.not.query + neither
    
    ## p-value calculation
    p.values <- matrix(
      unlist(
        lapply( row.names( query.and.sets ), function( k ) {
          .fisher_p(query.and.sets[k,], sets.not.query[k,], query.all[k,], 
                    query.colsum[k,]) 
        })), ncol=ncol(query), byrow=TRUE,
      dimnames=list(colnames(sets), colnames(query))
    )
    
    lor <- log((query.and.sets * neither) / (query.not.sets * sets.not.query))
    lor[query.not.sets == 0] <- Inf
    lor[sets.not.query == 0] <- Inf
    lor[query.and.sets == 0] <- 0
    
    ## store results
    results <- lapply( seq( ncol( query) ), function( g ) {
      res <- data.frame(
          set = colnames(sets),
          trend = ifelse(lor[,g] <= 0, "under", "over"),
          pval = p.values[,g ],
          padj = p.adjust( p.values[,g ], method="BH"),
          effect = round(.zScores( p.values[,g], direction=lor[,g]), 2),
          LOR = lor[,g],
          nSet = Matrix::colSums(abs(sets)),
          nFound = query.and.sets[,g ],
          signed = FALSE)
    })
    ## create separate result object for each query column
    if( length( results ) == 1 ){
      return( results[[ 1 ]])
    } else {
      return( results )
    }
}

.fisher_p <- function( x, y, m, n, relErr = 1 + 1e-7 ) {
  ## 'x' and 'y' are entries in the top two cells; 
  ## 'm' and 'n' are column totals.
  ## Code is excerpted from fisher.test, for efficiency. Note that 'support'
  ## varies in length with the input variables, so vectorization is only 
  ## possible via an mapply.
  mapply(
    function( x, y, m, n ) {
      k <- x + y
      lo <- max( 0, k - n )
      hi <- min( k, m )
      support <- lo:hi
      d <- dhyper( support, m, n, k, log = TRUE )
      d <- exp( d - max( d ) )
      d <- d / sum( d )
      sum( d[ d <= d[ x - lo + 1 ] * relErr ] )
    },
    x, y, m, n
  )
}

.zScores <- function(pval, direction=NULL, tails=2,
                     limit=.Machine$double.xmin) {
  if( !is.null( limit ) ){
    pval[which(pval < limit )] <- limit 
    ## set lower limit to avoid Inf/-Inf zscores
  }
  if( tails == 2){
    z <- qnorm( pval/2, lower.tail=FALSE )
  } else if( tails == 1){
    z <- qnorm(pval, lower.tail = FALSE)
  } else {
    stop( "Parameter 'tails' must be set to either 1 or 2.")
  }
  if ( !is.null( direction) ) {
    z <-  z * sign( direction )
  }
  z
}
