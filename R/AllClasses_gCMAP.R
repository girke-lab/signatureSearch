# Helper function authored by Brian Ripley, reproduced from
# http://r.789695.n4.nabble.com/checking-if-a-package-is-installed-td2340534.html
.f_checkpackage <- function(pkg) {
  tryCatch(require(pkg, character.only = TRUE), 
           error = function(e) FALSE) 
}

setClass("CMAPCollection",
         contains = "eSet",
         prototype = prototype(
           new("VersionedBiobase",
               versions=c(classVersion("eSet"), CMAPCollection="1.0.0")),
           phenoData = new("AnnotatedDataFrame",
             data=data.frame(),
             varMetadata=data.frame(
               labelDescription=character(0),
               channelDescription=factor())))
         )

setClass(
         "SignedGeneSet",
         contains = "GeneColorSet",
         validity = function( object ) {
           if (! all(levels( object@geneColor ) %in% c( "down", "up" ) ) )
             return( "Levels for 'geneColor'/'geneSign' factor must be 'down' and/or 'up'." )
           if ( any( is.na( object@geneColor ) ) )
             return( "All 'geneColor'/'geneSign' entries must be 'down' or 'up'." )
           return( TRUE )
         }
         )