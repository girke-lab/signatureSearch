## Changes in version 1.0.3 (2020-02-07)
+ Supported converting `feaResult` object to `enrichResult` object in the
`clusterProfiler` package so that the plotting functionalities in the latter 
package such as dotplots and gene-concept networks could be applied to the
FEA enrichment results
+ Supported subsetting specific columns (treatments) in the refdb (e.g. lincs) 
for GESS methods
+ Updated comp_fea_res function to reduce number of characters in description

## Changes in version 1.0.2 (2020-01-21)
+ Fix bug: the enrichment results from DSEA methods and some of TSEA methods
were added an aditional 'ont' column where the GO itmes were subsetted to the 
selected ontology

## Changes in version 1.0.1 (2019-11-10)
+ Support windows by not depending `gCMAP` package
+ Deal with HDF5 files with functions in `HDF5Array` package

## Changes in version 1.0.0 (2019-10-23)
+ Initial version 

## Changes in version 0.99.20 (2019-10-22)
+ Submitted to Bioconductor
+ Major changes
  - used HDF5 file to read and write the matrix in batches
  - data stored in ExperimentHub
