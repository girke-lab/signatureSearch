##' @import gCMAP
##' @import Biobase
##' @importFrom parallel mclapply
setMethod(
  "connectivity_score_raw",
  signature(experiment = "eSet", query = "CMAPCollection" ),
  function( experiment, query, element="z", keep.scores=FALSE) {
    if ( !( element %in% assayDataElementNames( experiment ) ) )
      stop( "Requested element name not found in data." )
    
    if( any( signed(query) == FALSE )) {
      stop("CMAPCollection contains unsigned GeneSets, which lack the information
about up-/down-regulated categories required to compute the connectivity score.")
    }

    ## rank data matrix in descending order
    data.matrix <- as( assayDataElement(experiment, element), "matrix" )

    ## subset objects to shared genes
    matched.features <- match( featureNames( experiment ), featureNames(query))
    matched.sets <- query[na.omit(matched.features),]
    
    ## extract scores for each gene set
    sets.up     <- mclapply( seq(ncol(matched.sets)),
                            function( x ) which(members( matched.sets )[ ,x ] == 1 ))
    
    sets.down     <- mclapply( seq(ncol(matched.sets)),
                              function( x ) which(members( matched.sets )[ ,x ] == -1 ))
    
    ## transform experiment to (reverse) ranks
    rank.matrix <- apply(data.matrix, 2, function(x) { length(x) - rank(x) + 1 } )
    
    ## calculate connectivity score
    raw.score <- apply( rank.matrix, 2, function( r ) {
      vapply(seq_along( sets.up ), function( n ) {
        .s( r[sets.up[[n]]], r[sets.down[[n]]], length( r ) )
      }, FUN.VALUE = numeric(1))
    })
    raw.score <- matrix(raw.score, ncol=ncol( experiment ))

    ## scale score across all tested query sets
    ## ThG: modification on next two lines:
    ## score <- apply(raw.score, 2, .S)
    score <- raw.score
    score <- matrix(score, ncol=ncol( experiment ))

    ## store raw per-gene expression scores
    if( keep.scores == TRUE) {
      gene.scores <- featureScores( experiment, query, element=element )
    } else { 
      gene.scores <- NA
             }
    ## store results
    results <- mclapply( seq( ncol( experiment ) ), function( x ) { ## x = data column
      
      if( ! all(is.na( gene.scores) )) {
        geneScores <- I( gene.scores[[x]])
      } else {
        geneScores <- NA
      }
  
      res <- data.frame(   set = sampleNames(query), 
                           trend = ifelse(score[,x] >=0, "up", "down"),
                           effect = score[,x],
                           nSet = Matrix::colSums( abs( members (query) ) ),
                           nFound = Matrix::colSums( abs( members (matched.sets) ) ),
                           pData(query))
     res
      })
    
    names( results ) <- sampleNames(experiment)
    
    ## return single CMAPResults of list of CMAPResults objects
    if( length( results ) == 1) {
      return ( results[[1]] )
    } else {
      return ( results )
    }
  }
          )

setMethod(
          "connectivity_score_raw",
          signature( experiment = "matrix",query = "CMAPCollection" ),
          function( experiment,query, ... ) {
            connectivity_score_raw( ExpressionSet(experiment), query, element="exprs" )
          }
          )

setMethod(
          "connectivity_score_raw",
          signature( experiment = "matrix" , query = "SignedGeneSet"),
          function( experiment, query, ...) {
            connectivity_score_raw( ExpressionSet(experiment), as(query, "CMAPCollection"), element="exprs")
          }
          )

setMethod(
          "connectivity_score_raw",
          signature( experiment = "eSet" ,query = "SignedGeneSet"),
          function( experiment, query, ...) {
            connectivity_score_raw( experiment, as(query, "CMAPCollection") )
          }
          )

setMethod(
  "connectivity_score_raw",
  signature( experiment = "matrix" , query = "GeneSetCollection"),
  function( experiment, query, ... ) {
    connectivity_score_raw( ExpressionSet(experiment), as(query, "CMAPCollection"), element="exprs")
  }
  )

setMethod(
  "connectivity_score_raw",
  signature( experiment = "eSet", query = "GeneSetCollection" ),
  function( experiment, query, ... ) {
    connectivity_score_raw(experiment, as(query, "CMAPCollection"), ...)
  }
  )

setMethod(
  "connectivity_score_raw",
  signature( experiment = "ANY", query = "GeneSet" ),
  function( experiment, query, ...) {
    stop("Connectivity score calculation requires gene sign information (up- / down- regulated gene categories).\n")
  }
  )


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
