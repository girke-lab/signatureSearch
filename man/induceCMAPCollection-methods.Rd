\name{induceCMAPCollection-methods}
\docType{methods}

\alias{induceCMAPCollection-methods}
\alias{induceCMAPCollection}

\alias{induceCMAPCollection,eSet-method}
\alias{induceCMAPCollection,matrix-method}
\title{  Methods for Function \code{induceCMAPCollection} in Package \pkg{gCMAP} }

\description{
This method defines a \code{\linkS4class{CMAPCollection}} by applying thresholds to an
element of an eSet-derived object. For example, applying
'induceCMAPCollection' to a matrix of z-scores stored in an NChannelSet,
gene sets can be defined for each of the sample columns stored in the
object. (See example section).
}

\usage{
\S4method{induceCMAPCollection}{eSet}(eset,element,lower=NULL,higher=NULL,sign.sets=TRUE)
\S4method{induceCMAPCollection}{matrix}(eset,element,...)
}

\arguments{
  \item{eset}{ An object derived from class \code{\linkS4class{eSet}},
    e.g. an NChannelSet}
  \item{element}{ A character string corresponding to the assayDataElementName of the 'eset' object to which the thresholds should be applied.}
  \item{lower}{ The lower threshold. If not 'NULL', genes with a score
    smaller than 'lower' will be included in the gene set with sign
  -1. At least one of 'lower' and 'higher' must be specified.}
  \item{higher}{ The 'higher' threshold. If not 'NULL', genes with a score
    larger than 'higher' will be included in the gene set with sign
    +1. At least one of 'lower' and 'higher' must be specified.}
  \item{sign.sets}{ Logical, indicating whether the 'signed' slot of the generated
    CMAPCollection should be set to 'TRUE' or 'FALSE'. This parameter
    should be set to 'FALSE' when the 'element' does not contain
    information about directionality, e.g. if it is a p-value.}
  \item{\dots}{ Any of the additional arguments detailed above.}
}

\value{
  A CMAPCollection object 
}

\keyword{methods}
\seealso{\code{\link{CMAPCollection}}
}
\examples{
\dontrun{
data(gCMAPData)
assayDataElementNames(gCMAPData)
cmap <- induceCMAPCollection(gCMAPData, element="z", lower=-2, higher=2)
cmap
notes(cmap)
}
}
