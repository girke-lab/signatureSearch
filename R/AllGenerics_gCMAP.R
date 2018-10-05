setGeneric("CMAPCollection",
           function(assayData,
                    phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                    featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                    annotation=character(),
                    protocolData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                    ...)
           standardGeneric("CMAPCollection"),
           signature="assayData")


setGeneric("induceCMAPCollection",
           def = function( eset, ... ) standardGeneric( "induceCMAPCollection" )
           )

setGeneric("setSizes",
           def = function( object) standardGeneric( "setSizes" )
)

setGeneric("members",
           def = function( object) standardGeneric( "members" )
           )

setGeneric(
           "signed",
           def = function( object) standardGeneric( "signed" )
           )

setGeneric("signed<-",
           def=function(x, value) standardGeneric("signed<-")
           )


setGeneric("mergeCollections",
           def = function( x, y) standardGeneric( "mergeCollections" )
           )

setGeneric("plot",
           def=function(x, y, ...) standardGeneric("plot")
           )

setGeneric("SignedGeneSet",
           def = function( type, ... ) standardGeneric( "SignedGeneSet" )
           )

setGeneric("upIds",
           def = function( object, ... ) standardGeneric( "upIds" )
           )

setGeneric("downIds",
           def = function( object, ... ) standardGeneric( "downIds" )
           )

setGeneric("geneSign",
           def = function( obj ) standardGeneric( "geneSign" )
           )

setGeneric("geneSign<-",
           def = function( object, value ) standardGeneric( "geneSign<-" )
           )

setGeneric("connectivity_score",
           def = function( experiment, query, ... ) standardGeneric( "connectivity_score" )
           )

setGeneric("fisher_score",
           def = function( query, sets, universe, ... ) standardGeneric( "fisher_score" )
           )

setGeneric("geneIndex",
           def = function( gene.sets, gene.ids, ...) standardGeneric( "geneIndex" )
           )

setGeneric("minSetSize",
           def = function( sets, ... ) standardGeneric( "minSetSize" )
           )

