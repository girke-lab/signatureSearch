\name{SignedGeneSet-class}

\Rdversion{1.1}

\docType{class}

\alias{SignedGeneSet-class}
\alias{downIds}
\alias{downIds,SignedGeneSet-method}
\alias{geneSign}
\alias{geneSign,SignedGeneSet-method}
\alias{geneSign<-}
\alias{geneSign<-,SignedGeneSet,character-method}
\alias{geneSign<-,SignedGeneSet,factor-method}
\alias{initialize,SignedGeneSet-method}
\alias{show,SignedGeneSet-method}
\alias{upIds}
\alias{upIds,SignedGeneSet-method}
\alias{mapIdentifiers,SignedGeneSet,GeneIdentifierType,GeneIdentifierType-method}
\alias{mapIdentifiers,SignedGeneSet,GeneIdentifierType,missing-method}
\alias{mapIdentifiers,SignedGeneSet,missing,ANY-method}
\alias{incidence,SignedGeneSet-method}
\alias{incidence,GeneSetCollection-method}

\title{Class \code{"SignedGeneSet"}}

\description{
  A simple extension of \code{\linkS4class{GeneColorSet}} which forces
  \code{geneColor} to be either \code{"down"} or \code{"up"} and which
  ignores \code{phenotype} and \code{phenotypeColor} slots.
}

\section{Objects from the Class}{

  Construct a \code{SignedGeneSet} with the \code{SignedGeneSet}
  constructor method, or with a call to \code{\link{new}}. Although
  \code{SignedGeneSet} derives from the more abstract
  \code{GeneColorSet}, not \code{phenotype} argument is required; if
  \code{phenotype} is supplied (or is present in a template object), it
  will be ignored.
  
}

\section{Slots}{
  See \code{\linkS4class{GeneColorSet}}. No additional slots are added.
}

\section{Extends}{
  Class \code{"\linkS4class{GeneColorSet}"}, directly.
  Class \code{"\linkS4class{GeneSet}"}, by class "GeneColorSet", distance 2.
}

\section{Methods}{

  Methods specific to \code{SignedGeneSet}:
  
  \describe{

    \item{downIds}{ \code{signature(object = "SignedGeneSet")}: retrieve
      \code{geneIds} entires for which \code{geneSign == "down"}. }

    \item{geneSign}{ \code{signature(obj = "SignedGeneSet")}: alias for
      \code{geneColor} slot. }

    \item{geneSign<-}{\code{signature(object = "SignedGeneSet", value = "character")}: alias for
      \code{geneColor} slot, converting to factor automatically. }

    \item{geneSign<-}{\code{signature(object = "SignedGeneSet", value = "factor")}: alias for
      \code{geneColor} slot. }

    \item{initialize}{\code{signature(.Object = "SignedGeneSet")}: on
      construction, checks for appropriate
      \code{geneSign}/\code{geneColor} values and sets \code{phenotype}
      and \code{phenotypeColor} to empty strings, since these are
      ignored. If no \code{geneSign}/\code{geneColor} values are
      supplied, \code{"up"} will be used by default. }

    \item{show}{\code{signature(object = "SignedGeneSet")}: same as for
      \code{\linkS4class{GeneColorSet}} but suppresses display of unused
      \code{phenotype} and \code{phenotypeColor} slots. }

    \item{upIds}{\code{signature(object = "SignedGeneSet")}: retrieve
      \code{geneIds} entires for which \code{geneSign == "up"}. }
    
    \item{mapIdentifiers}{{\code{signature(object = "SignedGeneSet")}}:
      Extends the 'mapIdentifiers' method implemented for GeneSets in
      the GSEABase package, but rejects target gene ids when multiple
      different (probe) identifiers with different gene signs (up /
      down) map to the same target.}
    
    \item{incidence}{{\code{signature(object = "SignedGeneSet")}}
      and
    }
    \item{incidence}{{\code{signature(object = "GeneSetCollection")}}:
      Mirror the 'incidence' method implemented for GeneSets in
      the GSEABase package, but returns sparseMatrix objects containing
      -1 / +1 to indicate up- and down-regulated gene members.}
  }

}

\keyword{classes}
\examples{
\dontrun{
gene.ids <- letters[1:10]
gene.signs <- rep(c("up","down"), each=5)
SignedGeneSet(gene.ids, geneSign=gene.signs, setName="set1")
}
}
