\name{CMAPCollection-class}
\Rdversion{1.1}
\docType{class}
\alias{CMAPCollection-class}

\alias{CMAPCollection}
\alias{CMAPCollection,environment-method}
\alias{CMAPCollection,Matrix-method}
\alias{CMAPCollection,matrix-method}
\alias{CMAPCollection,missing-method}

\alias{initialize,CMAPCollection-method}
\alias{geneIds,CMAPCollection-method}
\alias{setSizes,CMAPCollection-method}
\alias{members,CMAPCollection-method}
\alias{signed,CMAPCollection-method}
\alias{signed<-,CMAPCollection-method}
\alias{incidence,CMAPCollection-method}
\alias{upIds,CMAPCollection-method}
\alias{downIds,CMAPCollection-method}
\alias{incidence,list-method}
\alias{mergeCollections,CMAPCollection,CMAPCollection-method}

\alias{setSizes}
\alias{members}
\alias{signed}
\alias{signed<-}
\alias{mergeCollections}

\title{Class \code{"CMAPCollection"}}
\description{
  An extension of the eSet class for the efficient storage of (large) gene set collections.
}

\section{Objects from Class CMAPCollection}{
  
Objects can be created by calls of the form
\code{new("CMAPCollection", assayData, phenoData, featureData,
experimentData, annotation, protocolData, ...)}.  Alternatively, the
user-friendly 'CMAPCollection' method is available. 

The induceCMAPCollection function can be used to apply thresholds to
numerical scores stored in eSet-like objects and returns a CMAPCollecion
(see examples).

The CMAPCollection class is derived from the virtual
\code{\linkS4class{eSet}} class. The assayData slot contains information
about the membership of genes (rows) in gene sets (columns) in the form
of an incidence matrix. The incidence matrix, accessible through the
'members' method, is a 'sparseMatrix' object, in which 1 / -1 entries
identify gene set membership of up- and downregulated genes,
respecively.

As opposed to the well-established GeneSetCollection class defined in
the GSEABase package, the CMAPCollection class stores gene set
membership in a matrix format, allowing direct access to individual gene
sets as well as the relationships between different sets. The incidence
matrix offers memory efficient storage of large gene set collection and
can directly be used in matrix-based gene set analyses.

Through direct extension of the virtual eSet class, featureData and phenoData
slots are available for storage of gene- and gene-set annotation,
respectively. The column 'signed' in the phenoData slot indicates
whether the different gene sets ( columns ) should be considered to be
signed to disambiguate cases in which all gene set members are
identified by a +1 entry. In this case, 'signed' = TRUE indicates that these genes should be
considered upregulated members of the set (and no down-regulated members
were identified / stored). If 'signed' = FALSE, no information about
directionality is available, e.g. gene set members can be either up- or
downregulated.

}
\section{Slots}{
  \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"}  }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"}  }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{geneIds}{\code{signature(object = "CMAPCollection")}:
      Returns a list of gene identifiers, with one list entry for each
      column of the assayDataSlot 'members'. }
    \item{members}{\code{signature(object = "CMAPCollection")}:
      Returns the number of gene members in each gene set stored in the collection.
      For signed sets, also the number of up-/down-regulated members is returned.}
    \item{members}{\code{signature(object = "CMAPCollection")}:
      Returns the coincidence matrix as stored in the assayData slot of the
      CMAPCollection as a sparseMatrix object (rows=genes, columns=gene sets). }
    \item{signed}{\code{signature(object = "CMAPCollection")}:
      Returns the 'signed' column of the phenoData slot, indicating
      whether gene sets should be considered signed (TRUE) or un-signed (FALSE).}
    \item{signed<-}{\code{signature(x = "CMAPCollection")}:
      Replacement method for the 'signed' column of phenoData.}
    \item{minSetSize}{\code{signature(sets =
	"CMAPCollection")}:
      Filter CMAPcollection for minimum number of set members.}
    \item{incidence}{\code{signature(x = "CMAPCollection")}:
      Returns in the transpose of the coincidence matrix stored in the
      assayData slot, mirroring the definition of 'incidence' for
      GeneSetCollections as defined in the GSEABase package.}
    \item{mergeCollections}{\code{signature(x = "CMAPCollection", y =
	"CMAPCollection")}: Combines two CMAPCollections into one.}
    \item{upIds}{\code{signature(x = "CMAPCollection")}: Returns the gene identifiers of all up-regulated gene set members (sign = 1).}
    \item{downIds}{\code{signature(x = "CMAPCollection")}: Returns the gene identifiers of all down-regulated gene set members (sign = -1).}
  }
  }

\author{
Thomas Sandmann, sandmann.thomas@gene.com
}
\note{
The CMAPCollections class supports coercion from / to GeneSet and
GeneSetCollection objects defined by the GSEABase package, as well as
the SignedGeneSet derivative introduced by the gCMAP package itself.
}

\seealso{
  \code{\link{induceCMAPCollection}},
  \code{\linkS4class{GeneSetCollection}}, \code{\linkS4class{SignedGeneSet}}
}
\examples{
\dontrun{
## empty CMAPCollection
new("CMAPCollection")

## CMAPCollection from matrix
mat <- matrix( sample( c(-1,0,1), 100, replace=TRUE), ncol=10)
cmap <- CMAPCollection( mat )
members( cmap )

## CMAPCollection induced from NChannelSet
data( gCMAPData )
assayDataElementNames( gCMAPData )

cmap <- induceCMAPCollection(gCMAPData, "z", lower=-2, higher=2)
cmap
setSizes( cmap )
pData(cmap)
signed(cmap) <- c(TRUE, FALSE, TRUE)
signed(cmap)
head(members(cmap))

out <- fisher_score(cmap[,1], cmap, universe = featureNames( cmap))
out
}
}

\keyword{classes}
