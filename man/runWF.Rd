% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runWF.R
\name{runWF}
\alias{runWF}
\title{Run the Entire GESS/FEA Workflow}
\usage{
runWF(
  Signature,
  cellInfo,
  PertColName = "pert_iname",
  drug,
  refdb,
  gess_method = "LINCS",
  fea_method = "dup_hyperG",
  N_gess_drugs = 150,
  env_dir = ".",
  tau = FALSE,
  Nup = 150,
  Ndown = 150,
  higher = 1,
  lower = -1,
  method = "spearman",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 5,
  maxGSSize = 500,
  runFEA = TRUE,
  GenerateReport = TRUE
)
}
\arguments{
\item{Signature}{The signature to perform signatureSearching on.}

\item{cellInfo}{A data table containing information about the cell types 
that make up the signature search database. This table must contain a column
labeled "cell_id" that contains matching cell type names in the signature 
search database.  Details about this file  can be found 
in the \code{cell_info2} table after load the `signatureSearch` package
and running `data("cell_info2")`}

\item{PertColName}{The name of the column from the signature search results 
table used to summarize results by.}

\item{drug}{character(1) representing query drug name (e.g. vorinostat). 
This query drug should be included in the \code{refdb}}

\item{refdb}{character(1), one of "lincs", "lincs_expr", "cmap", "cmap_expr",
or path to the HDF5 file built from \code{\link{build_custom_db}} function}

\item{gess_method}{character(1), one of "LINCS", "CORsub", "CORall", 
"Fisher", "CMAP", "gCMAP". When \code{gess_method} is "CORsub" or "CORall",
only "lincs_expr" or "cmap_expr" databases are supported.}

\item{fea_method}{character(1), one of "dup_hyperG", "mGSEA", "mabs", 
"hyperG", "GSEA"}

\item{N_gess_drugs}{number of unique drugs in GESS result used as input of
FEA}

\item{env_dir}{character(1), directory under which the result environment 
located. The default is users current working directory in R session, can
be checked via \code{getwd()} command in R}

\item{tau}{TRUE or FALSE indicating whether to compute Tau scores if 
\code{gess_method} is set as 'LINCS'}

\item{Nup}{integer(1). Number of most up-regulated genes to be subsetted 
for GESS query when \code{gess_method} is CMAP, LINCS or CORsub}

\item{Ndown}{integer(1). Number of most down-regulated genes to be subsetted 
for GESS query when \code{gess_method} is CMAP, LINCS or CORsub}

\item{higher}{numeric(1), it is defined when gess_method argument is 'gCMAP'
or 'Fisher' representing the 'upper' threshold of subsetting genes with a 
score larger than 'higher'}

\item{lower}{numeric(1), it is defined when gess_method argument is 'gCMAP'
or 'Fisher' representing the 'lower' threshold of subsetting genes}

\item{method}{One of 'spearman' (default), 'kendall', or 'pearson', 
indicating which correlation coefficient to use}

\item{pvalueCutoff}{double, p-value cutoff for FEA result}

\item{qvalueCutoff}{double, qvalue cutoff for FEA result}

\item{minGSSize}{integer, minimum size of each gene set in annotation system}

\item{maxGSSize}{integer, maximum size of each gene set in annotation system}

\item{runFEA}{Logical value indicating if FEA analysis is performed.}

\item{GenerateReport}{Logical value indicating if a report is generated.}
}
\value{
list object containing GESS/FEA result tables
}
\description{
This function runs the entire GESS/FEA workflow when providing 
the query drug and cell type, as well as selecting the reference database 
(e.g. 'cmap' or 'lincs'), defining the specific GESS and FEA methods. 
In this case, the query GES is drawn from the reference database. 
The N (defined by the `N_gess_drugs` argument) top ranking hits in the GESS 
tables were then used for FEA where three different annotation systems were 
used: GO Molecular Function (GO MF), GO Biological Process (GO BP) and 
KEGG pathways. 

The GESS/FEA results will be stored in a list object in R session. 
A working environment named by the use case will be created under users 
current working directory or under other directory defined by users.
This environment contains a \code{results} folder where the GESS/FEA 
result tables were written to. The working environment also contains
a template Rmd vignette as well as a rended HTML report, users could make
modifications on the Rmd vignette as they need and re-render it to generate 
their HTML report.
}
\examples{
library(signatureSearch)
library(ExperimentHub); library(rhdf5)
eh <- ExperimentHub()
cmap <- eh[["EH3223"]]; cmap_expr <- eh[["EH3224"]]
lincs <- eh[["EH3226"]]; lincs_expr <- eh[["EH3227"]]
lincs2 <- eh[["EH7297"]]
h5ls(lincs2)
db_path <- system.file("extdata", "sample_db.h5", package = "signatureSearch")
library(SummarizedExperiment);
library(HDF5Array)
sample_db <- SummarizedExperiment(HDF5Array(db_path, name="assay"))
rownames(sample_db) <- HDF5Array(db_path, name="rownames")
colnames(sample_db) <- HDF5Array(db_path, name="colnames")
query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
query <- as.numeric(query_mat); names(query) <- rownames(query_mat)
upset <- head(names(query[order(-query)]), 150)
head(upset)
downset <- tail(names(query[order(-query)]), 150)
head(downset)
runWF(Signature = list(upset=upset, downset=downset),
      cellInfo = cellInfo2,
      PertColName = "pert_iname",
      drug = "vorinostat",
      refdb = lincs2,
      gess_method="LINCS",
      fea_method="dup_hyperG",
      N_gess_drugs=150,
      env_dir="./GESSWFResults",
      tau=FALSE,
      Nup=150,
      Ndown=150,
      higher=1,
      lower=-1,
      method="spearman",
      pvalueCutoff=1,
      qvalueCutoff=1,
      minGSSize=5,
      maxGSSize=500,
      runFEA=TRUE,
      GenerateReport= TRUE)
}
