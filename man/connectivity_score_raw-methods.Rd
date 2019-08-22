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

\description{ A modified version of the original \code{connectivity_score} 
method in `gCMAP` package, which stores the raw connectivity scores before 
scaling in the result.}
\details{
The reference database for \code{gess_gcmap} method is an hdf5 file storing
gene expression matrix from a large number of compound treatments. The matrix
in hdf5 file is read into R by chunks and induced to a 
\code{\link[gCMAP]{CMAPCollection}} object by 
\code{\link[gCMAP]{induceCMAPCollection}} method and setting higher and lower
cutoffs. The query matrix in the 'qSig' argument of the \code{gess_gcmap}
function and the induced \code{CMAPCollection} object are used as the input
for 'connectivity_score_raw' function to compute raw connectivity scores. 
The scaling is applied after all the signatures in the reference database are
searched.
}

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
  
}

\seealso{
\code{\link[gCMAP]{connectivity_score}}; 
\code{\link[gCMAP]{CMAPCollection}}; 
\code{\link[gCMAP]{induceCMAPCollection}}
}

