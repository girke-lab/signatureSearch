\name{connectivity_score}

\alias{ks}
\alias{s}
\alias{S}

\alias{connectivity_score}
\alias{connectivity_score,eSet,CMAPCollection-method}
\alias{connectivity_score,matrix,CMAPCollection-method}
\alias{connectivity_score,matrix,SignedGeneSet-method}
\alias{connectivity_score,eSet,SignedGeneSet-method}
\alias{connectivity_score,eSet,GeneSetCollection-method}
\alias{connectivity_score,matrix,GeneSetCollection-method}
\alias{connectivity_score,ANY,GeneSet-method}

\title{Broad CMAP gene set enrichment metrics}

\description{ A method for computing Broad CMAP connectivity scores, as
  described in the reference below. Supporting functions used for
  computation are also described. }

\usage{
\S4method{connectivity_score}{eSet,CMAPCollection}(experiment, query, element="z", keep.scores=FALSE)

\S4method{connectivity_score}{matrix,CMAPCollection}(experiment, query, ...)

\S4method{connectivity_score}{eSet,SignedGeneSet}(experiment, query, ...)

\S4method{connectivity_score}{matrix,SignedGeneSet}(experiment, query, ...)

\S4method{connectivity_score}{eSet,GeneSetCollection}(experiment, query,
...)

\S4method{connectivity_score}{matrix,GeneSetCollection}(experiment, query,...)

\S4method{connectivity_score}{ANY,GeneSet}(experiment, query, ...)
}

\arguments{

  \item{experiment}{ An \code{\linkS4class{eSet}} or matrix object to query.}
  
  \item{query}{ A \code{\linkS4class{CMAPCollection}}, \code{\linkS4class{SignedGeneSet}}, or
    \code{\linkS4class{GeneSetCollection}} object containing signed gene
    sets with which to query the \code{experiment} object. }

  \item{element}{ Character string specifying which element of a
    multi-channel eSet to access for determining tag rank? }
  
\item{keep.scores}{ Logical: keep gene-level scores for all gene sets (Default: TRUE) ? The size of the generated CMAPResults object increases with the number of contained gene sets. For very large collections, consider setting this parameter to 'FALSE' to conserve memory.}

\item{\dots}{ Additional arguments passed on to downstream functions.}
}

\value{

  \item{\code{connectivity_score}}{For the
    \code{\linkS4class{SignedGeneSet}} method, a vector of \code{s}
    scores, one per instance in \code{experiment}. For the
    \code{\linkS4class{GeneSetCollection}} method, a matrix, with one
    row per instance in \code{experiment} and one column per query set. }

  \item{\code{ks}}{ A signed Kolmogorov-Smirnov type statistic based on
    the position of the ranks \code{V} in the vector \code{1:n}. }

  \item{\code{s}}{ A difference of \code{ks} values for \code{V_up}
    vs. \code{V_down}, or 0 if both yield the same sign. }

  \item{\code{S}}{ A vector of signed, rescaled scores. After rescaling,
    1 corresponds to the maximum positive \code{s} score, and -1, to the
    minimum negative \code{s} score. \code{S} is typically used to
    produce the red-grey-green instance heat maps from the reference
    below. }
  
}

\note{Note that as defined by Lamb et al., \code{ks} is not
  symmetric. For \code{n = 100}, for example, \code{ks(1,100)} is .99
  while \code{ks(100,100} is -1. A further consequence of the Lamb et
  al. definitions is that the maximum possible score for a perfect
  positive match depends on query set size. See the example
  below. Although these properties are not desirable, the intention here
  is to exactly reproduce the Lamb et al. statistic.}

\references{ Lamb, J. et al. (2006). The Connectivity Map: using
gene-expression signatures to connect small molecules, genes, and
disease. \emph{Science} 313:1929. Notation for \code{ks}, \code{s}, and
\code{S} closely follows the Supporting Online Material there. }

\examples{
\dontrun{
data(gCMAPData)

  ## induce CMAPCollection from z-scores
  sets <- induceCMAPCollection(gCMAPData, "z", lower=-3, higher=3)

  ## Broad CMAP KS scoring: one z-score column
  connectivity_score(gCMAPData[,1], sets, element="z")

  ## multiple z-score columns, results are returned as a list
  connectivity_score(gCMAPData, sets)
}
}
