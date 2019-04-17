
# _signatureSearch_: Discovering novel modes of action of bioactive compounds

# Introduction

This project is about optimizing signature search and enrichment methods for 
the discovery of novel modes of action (MOA) of bioactive compounds from 
reference databases, such as LINCS, containing the genome-wide gene expression 
signatures (GESs) from tens of thousands of drug and genetic perturbations. 
The methods can be divided into two major classes. First, gene expression 
signature search (GESS) methods are used to identify drugs that induce GESs 
similar to those of query GESs of interest. The queries can be drug- or 
disease-related GESs. Since the MOA of most drugs in the corresponding 
reference databases are known, the resulting associations are useful to gain 
insights into pharmacological and disease mechanisms and to develop novel drug 
repurposing approaches. Second, functional enrichment analysis (FEA) methods 
using Gene Ontology (GO) or pathway annotations have been developed to 
functionally interpret the vast number of GESS results. The latter are 
composed of lists of drugs ranked by the similarity metric of the corresponding 
GESS method making the functional interpretation of their top ranking drugs 
challenging. Importantly, the FEA methods developed by this study also support 
the reconstruction of drug-target networks to guide the interpretation of the 
results.

# Installation and Loading

`signatureSearch` is a Bioconductor package and can be installed through 
`BiocManager::install()`

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("signatureSearch", version = "3.9")
```
Or it can be direcly installed from GitHub by
```r
devtools::install_github("yduan004/signatureSearch")
```

After the package is installed, it can be loaded into *R* workspace by
```r
library(signatureSearch)
```
For detailed description of the package, please refer to the vignette by
```r
browseVignettes("signatureSearch")
```


