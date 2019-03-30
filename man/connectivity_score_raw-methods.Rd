\name{connectivity_score_raw}

\alias{ks}
\alias{s}
\alias{S}

\alias{connectivity_score_raw}
\alias{connectivity_score_raw,eSet,CMAPCollection-method}
\alias{connectivity_score_raw,matrix,CMAPCollection-method}
\alias{connectivity_score_raw,matrix,SignedGeneSet-method}
\alias{connectivity_score_raw,eSet,SignedGeneSet-method}
\alias{connectivity_score_raw,eSet,GeneSetCollection-method}
\alias{connectivity_score_raw,matrix,GeneSetCollection-method}
\alias{connectivity_score_raw,ANY,GeneSet-method}

\title{Broad CMAP gene set enrichment metrics}

\description{ A modified version of the original connectivity_score method in 
             `gCMAP` package, which stores the raw connectivity scores 
             before scaling in the result.}

\usage{
\S4method{connectivity_score_raw}{eSet,CMAPCollection}(experiment, query, 
element="z", keep.scores=FALSE)

\S4method{connectivity_score_raw}{matrix,CMAPCollection}(experiment, query, ...)

\S4method{connectivity_score_raw}{eSet,SignedGeneSet}(experiment, query, ...)

\S4method{connectivity_score_raw}{matrix,SignedGeneSet}(experiment, query, ...)

\S4method{connectivity_score_raw}{eSet,GeneSetCollection}(experiment, query,
...)

\S4method{connectivity_score_raw}{matrix,GeneSetCollection}(experiment, 
query,...)

\S4method{connectivity_score_raw}{ANY,GeneSet}(experiment, query, ...)
}

\arguments{

  \item{experiment}{ An \code{eSet} or matrix object to query.}
  
  \item{query}{ A \code{CMAPCollection}, \code{SignedGeneSet}, or
    \code{GeneSetCollection} object containing signed gene
    sets with which to query the \code{experiment} object. }

  \item{element}{ Character string specifying which element of a
    multi-channel eSet to access for determining tag rank? }
  
\item{keep.scores}{ Logical: keep gene-level scores for all gene sets 
(Default: TRUE) ? The size of the generated CMAPResults object increases with 
the number of contained gene sets. For very large collections, consider setting 
this parameter to 'FALSE' to conserve memory.}

\item{\dots}{ Additional arguments passed on to downstream functions.}
}

\value{

  \item{\code{connectivity_score_raw}}{For the
    \code{SignedGeneSet} method, a vector of \code{s}
    scores, one per instance in \code{experiment}. For the
    \code{GeneSetCollection} method, a matrix, with one
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

\seealso{\code{\link[gCMAP]{connectivity_score}}}
