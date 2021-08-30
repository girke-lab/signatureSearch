
# _signatureSearch_: Environment for Gene Expression Searching and Functional Enrichment Analysis
[![platforms](http://www.bioconductor.org/shields/availability/release/signatureSearch.svg)](http://www.bioconductor.org/packages/devel/bioc/html/signatureSearch.html#archives)
[![rank](http://www.bioconductor.org/shields/downloads/devel/signatureSearch.svg)](http://bioconductor.org/packages/stats/bioc/signatureSearch/)
[![posts](http://www.bioconductor.org/shields/posts/signatureSearch.svg)](https://support.bioconductor.org/t/signatureSearch/)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/signatureSearch.svg)](http://www.bioconductor.org/packages/devel/bioc/html/signatureSearch.html#since)
[![build](http://www.bioconductor.org/shields/build/devel/bioc/signatureSearch.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/signatureSearch/)
[![updated](http://www.bioconductor.org/shields/lastcommit/devel/bioc/signatureSearch.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/signatureSearch/)
[![dependencies](http://www.bioconductor.org/shields/dependencies/devel/signatureSearch.svg)](http://www.bioconductor.org/packages/devel/bioc/html/signatureSearch.html#since)

# Introduction

The _signatureSearch_ package implements algorithms and data structures for
performing gene expression signature (GES) searches, and subsequently
interpreting the results functionally with specialized enrichment methods.
These utilities are useful for studying the effects of genetic, chemical and
environmental perturbations on biological systems. Specifically, in drug
discovery they can be used for identifying novel modes of action (MOA) of
bioactive compounds from reference databases such as LINCS containing the
genome-wide GESs from tens of thousands of drug and genetic perturbations
(Subramanian et al. 2017). A typical GES search (GESS) workflow can be divided into
two major steps.  First, GESS methods are used to identify
perturbagens such as drugs that induce GESs similar to a query GES of interest.
The queries can be drug-, disease- or phenotype-related GESs. Since the
MOAs of most drugs in the corresponding reference databases are known, the
resulting associations are useful to gain insights into pharmacological and/or
disease mechanisms, and to develop novel drug repurposing approaches. Second,
specialized functional enrichment analysis (FEA) methods using annotations
systems, such as Gene Ontologies (GO), pathways or Disease Ontologies (DO),
have been developed and implemented in this package to efficiently interpret
GESS results. The latter are usually composed of lists of perturbagens (_e.g._
drugs) ranked by the similarity metric of the corresponding GESS method.
Finally, network resconstruction functionalities are integrated for visualizing
the final results, _e.g._ in form of drug-target networks. For each GESS and FEA
step, several alternative methods have been implemented in _signatureSearch_ to
allow users to choose the best possible workflow configuration for their
research application. 

# Vignette
The vignette of this package is available at [here](https://www.bioconductor.org/packages/release/bioc/vignettes/signatureSearch/inst/doc/signatureSearch.html)

# Installation and Loading

`signatureSearch` is a R/Bioconductor package and can be installed using 
`BiocManager::install()`.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("signatureSearch")
```
To obtain the most recent updates immediately, one can install it directly from 
GitHub as follows.
```r
devtools::install_github("yduan004/signatureSearch", build_vignettes=TRUE)
```

After the package is installed, it can be loaded into an R session as follows.
```r
library(signatureSearch)
```
For detailed description of the package, please refer to the vignette by running
```r
browseVignettes("signatureSearch")
```


