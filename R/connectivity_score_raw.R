##' @import gCMAP
##' @import Biobase
connectivity_score_raw <- function(experiment, query) {
    if(any(signed(query) == FALSE)) {
      stop("CMAPCollection contains unsigned GeneSets, which lack the
       information about up-/down-regulated categories required to
       compute the connectivity score.")
    }
    data.matrix <- experiment
    
    ## subset objects to shared genes
    matched.features <- match(rownames(experiment), rownames(query))
    matched.sets <- query[na.omit(matched.features),]
    
    ## extract scores for each gene set
    sets.up <- lapply(seq(ncol(matched.sets)),
                        function(x) which(members(matched.sets)[ ,x ] == 1))
    
    sets.down <- lapply(seq(ncol(matched.sets)),
                        function(x) which(members(matched.sets)[ ,x] == -1))
    
    ## transform experiment to (reverse) ranks
    rank.matrix <- apply(data.matrix, 2, function(x) {length(x)-rank(x)+1})
    
    ## calculate connectivity score
    raw.score <- apply(rank.matrix, 2, function(r) {
      vapply(seq_along(sets.up), function(n) {
        .s(r[sets.up[[n]]], r[sets.down[[n]]], length(r))
      }, FUN.VALUE = numeric(1))
    })
    raw.score <- matrix(raw.score, ncol=ncol(experiment))

    ## scale score across all tested query sets
    ## ThG: modification on next two lines:
    ## score <- apply(raw.score, 2, .S)
    score <- raw.score
    score <- matrix(score, ncol=ncol(experiment))
    ## store results
    results <- data.frame(set = sampleNames(query), 
                    trend = ifelse(score[,1] >=0, "up", "down"),
                    effect = score[,1],
                    nSet = colSums(as.matrix(abs(members(query)))),
                    nFound = colSums(as.matrix(abs(members(matched.sets)))),
                    pData(query))
}

.ks <- function( V, n ) {
  t <- length( V )
  
  if( t == 0 )  {
    return( 0 )
  } else {
    
    if ( is.unsorted( V ) )
      V <- sort( V )
    d <- seq_len(t) / t - V / n
    a <- max( d )
    b <- -min( d ) + 1 / t
    ifelse( a > b, a, -b )
  }
}

.s <- function( V_up, V_down, n ) {
  ks_up <- .ks( V_up, n )
  ks_down <- .ks( V_down, n )
  ifelse( sign( ks_up ) == sign( ks_down ), 0, ks_up - ks_down )
}

.S <- function( scores ) {
  p <- max( scores )
  q <- min( scores )
  ifelse(
         scores == 0,
         0,
         ifelse( scores > 0, scores / p, -scores / q )
         )
}
