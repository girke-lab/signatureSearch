
# _signatureSearch_: Discovering novel modes of action of bioactive compounds

# Background

This project is about optimizing signature search and enrichment methods for the discovery of novel modes of action (MOA) of bioactive compounds from reference databases, such as LINCS, containing the genome-wide gene expression signatures (GESs) from tens of thousands of drug and genetic perturbations. The methods can be divided into two major classes. First, gene expression signature search (GESS) methods are used to identify drugs that induce GESs similar to those of query GESs of interest. The queries can be drug- or disease-related GESs. Since the MOA of most drugs in the corresponding reference databases are known, the resulting associations are useful to gain insights into pharmacological and disease mechanisms and to develop novel drug repurposing approaches. Second, functional enrichment analysis (FEA) methods using Gene Ontology (GO) or pathway annotations have been developed to functionally interpret the vast number of GESS results. The latter are composed of lists of drugs ranked by the similarity metric of the corresponding GESS method making the functional interpretation of their top ranking drugs challenging. Importantly, the FEA methods developed by this study also support the reconstruction of drug-target networks to guide the interpretation of the results.

@Lamb2006-du generated a GES database, initially including 164 drugs screened against four mammalian cell lines [@Lamb2006-du]. A few years later it was extended to 1,309 drugs and eight cell lines. In 2017, the database was increased to 19,811 drugs by the Library of Network-Based Cellular Signatures (LINCS) Consortium [@Subramanian2017-fu], while the number of cell types represented in the database has increased to over 70 normal and cancer cell lines. The number of compound dosages and time points considered in the assays has also been increased by 10-20 fold. The initial database used Affymetrix Gene Chips as expression platform. To scale from a few thousand to many hundred thousand GESs, the LINCS Consortium uses now the more economic L1000 assay. This bead-based technology is a low cost, high-throughput reduced representation expression profiling assay. It measures the expression of 978 landmark genes and 80 control genes by detecting fluorescent intensity of beads after capturing the ligation-mediated amplification products of mRNAs [@Peck2006-rf]. The expression of 11,350 additional genes is imputed from the landmark genes by using as training data a collection of 12,063 Affymetrix gene chips [@Edgar2002-di]. The substantial scale-up of the LINCS project provides now many new opportunities to explore MOAs for a vast number of small molecules.

# Terminology

Gene Expression Signatures (GESs): can be gene expression profiles (GEPs), which can be genome-wide profile from differential expression analysis ($log_2$ ratios or z-scores) or gene expression intensity values from perturbagen treatment in cell culture, or differentially expressed gene sets (DEGs) from a treatment.

# Reference Database

The CMAP and LINCS differential expression databases can be built with `build_db` funciton or directly downloaded at [CMAP](http://biocluster.ucr.edu/~yduan004/CMAP_db/cmap.tar.gz), [LINCS](http://biocluster.ucr.edu/~yduan004/LINCS_db/lincs42.tar.gz). The databases are stored as HDF5 backed `SummarizedExperiment` object. Untar the downloaded file by `tar -xzvf file.tar.gz` command, then load the `SummarizedExperiment` object via `loadHDF5SummarizedExperiment` function. The 'assays' slot of the loaded SummarizedExperiment object represents the $log_2$ ratios or z-scores generated from differential expression (DE) analysis. 

The CMAP and LINCS expression databases can also be built with `build_db` funciton or directly downloaded at [CMAP_expr](http://biocluster.ucr.edu/~yduan004/CMAP_db/cmap_expr.tar.gz), [LINCS_expr](http://biocluster.ucr.edu/~yduan004/LINCS_db/lincs42_expr.tar.gz). They stores the mean experssion values of drug treatment samples in cells.

The custom database can also be built via `build_db` function if a `data.frame` representing genoime-wide GEPs of compound or genetic treatments in cells is provided.

# Methods for GESS

As an example of running GESS/FEA workflow, 95 GEPs randomly sampled from LINCS DE database and 5 GEPs from `HDAC inhibitors` (`vorinostat`, `rhamnetin`, `trichostatin A`, `valproic acid`, and `HC toxin`) treatments in MCF7 cell were combined to generate a sample reference database. The query signature is the vorinostat treatment drawn from the sample database.  

## Load database and get query signature
```r
library(signatureSearch)
db_dir <- system.file("extdata", "sample_db", package = "signatureSearch")
sample_db <- loadHDF5SummarizedExperiment(db_dir)
# get "vorinostat__MCF7__trt_cp" signature drawn from sample database
query_mat <- as.matrix(assay(sample_db[,"vorinostat__MCF7__trt_cp"]))
query = as.numeric(query_mat); names(query) = rownames(query_mat)
upset <- head(names(query[order(-query)]), 150)
downset <- tail(names(query[order(-query)]), 150)
```

## Searching with CMAP method for GESS

Generate "qSig" object for GESS_CMAP method and search against reference database
```r
qsig_cmap <- qSig(qsig = list(upset=upset, downset=downset), gess_method = "CMAP", refdb = sample_db)
cmap <- gess_cmap(qSig=qsig_cmap, chunk_size=5000)
cmap
result(cmap)
```

## Searching with LINCS method for GESS

Generate "qSig" object for GESS_LINCS method and search against reference database
```r
qsig_lincs <- qSig(qsig = list(upset=upset, downset=downset), gess_method = "LINCS", refdb = sample_db)
lincs <- gess_lincs(qsig_lincs, sortby="NCS")
result(lincs)
#saveRDS(lincs, "~/insync/project/lincs_against_sample_db.rds")
```

## Searching with gCMAP method for GESS

Generate "qSig" object for GESS_gCMAP method and search against reference database
```r
qsig_gcmap <- qSig(qsig = query_mat, gess_method = "gCMAP", refdb = sample_db)
gcmap <- gess_gcmap(qsig_gcmap, higher=1, lower=-1)
result(gcmap)
```

## Searching with Fisher method for GESS

Generate "qSig" object for GESS_fisher method and search against reference database
```r
qsig_fisher <- qSig(qsig = query_mat, gess_method = "Fisher", refdb = sample_db)
fisher <- gess_fisher(qSig=qsig_fisher, higher=1, lower=-1)
result(fisher)
```

## Searching with Spearman correlation method for GESS

Generate "qSig" object for GESS_Cor method and search against reference database

### Genome-wide Spearman correlation
```r
qsig_sp <- qSig(qsig = as.matrix(query), gess_method = "Cor", refdb = sample_db)
sp <- gess_cor(qSig=qsig_sp, method="spearman")
result(sp)
```

### Spearman_sub correlatioin
```r
# Subset z-scores of 150 up and down gene sets from "vorinostat__MCF7__trt_cp" signature.
query_mat_sub <- as.matrix(query_mat[c(upset, downset),])
qsig_spsub <- qSig(qsig = query_mat_sub, gess_method = "Cor", refdb = sample_db)
spsub <- gess_cor(qSig=qsig_spsub, method="spearman")
result(spsub)
```

# Methods for FEA

Choose GESS result from LINCS method as input for the downstream functional enrichment analysis.

## dup_hyperG method for TSEA

Subset top 10 ranking drugs in GESS result, get their target set (with duplication) as query for duplication support hypergeometric test for TSEA.
```r
drugs <- unique(result(lincs)$pert[1:10])
# GO annotation system
dup_hyperG_res <- tsea_dup_hyperG(drugs = drugs, universe = "Default", type = "GO", ont="MF", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff = 0.1, minGSSize = 10, maxGSSize = 500)
dup_hyperG_res
result(dup_hyperG_res)
# KEGG annotation system
dup_hyperG_k_res <- tsea_dup_hyperG(drugs = drugs, universe = "Default", type = "KEGG", pvalueCutoff=0.1, pAdjustMethod="BH", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500)
dup_hyperG_k_res
result(dup_hyperG_k_res)
```

## m_GSEA method for TSEA

Subset top 10 ranking drugs in GESS result as query, m_GSEA method internally get their target set and turn it to scored ranked target list as query for TSEA. 
The scores represent weight of targets in the target set.
```r
mgsea_res <- tsea_mGSEA(drugs=drugs, type="GO", ont="MF", exponent=1, nPerm=1000, pvalueCutoff=1, minGSSize=5)
mgsea_res
result(mgsea_res)
# KEGG annotation system
mgsea_k_res <- tsea_mGSEA(drugs=drugs, type="KEGG", exponent=1, nPerm=1000, pvalueCutoff=1, minGSSize=2)
mgsea_k_res
result(mgsea_k_res)
```

## mabs method for TSEA

Subset top 10 ranking drugs in GESS result as query, mabs method internally get their target set and ranked target list with scores, which represent weight of targets. 
```r
# GO annotation system
mabs_res <- tsea_mabs(drugs=drugs, type="GO", ont="MF", nPerm=1000, pvalueCutoff=0.05, minGSSize=5)
result(mabs_res)
# KEGG annotation system
mabs_k_res <- tsea_mabs(drugs=drugs, type="KEGG", nPerm=1000, pvalueCutoff=0.05, minGSSize=5)
result(mabs_k_res)
```

## hyperG method for DSEA

Subset top 10 ranking drugs in GESS result as query for hypergeometric test for DSEA.
```r
drugs <- unique(result(lincs)$pert[1:10])
# GO annotation system
hyperG_res <- dsea_hyperG(drugs = drugs, type = "GO", ont="MF", pvalueCutoff=0.05, qvalueCutoff = 0.1, minGSSize = 10)
hyperG_res
result(hyperG_res)
# KEGG annotation system
hyperG_k_res <- dsea_hyperG(drugs = drugs, type = "KEGG", pvalueCutoff=0.2, qvalueCutoff = 0.2, minGSSize = 10)
hyperG_k_res
result(hyperG_k_res)
```

## GSEA method for DSEA

Use ranked drug list in GESS result as GSEA input, the scores are similarity metrics of the corresponding GESS methods. Zeros are removed.
```r
dl <- abs(result(lincs)$NCS); names(dl) <- result(lincs)$pert
dl <- dl[dl>0]
dl <- dl[!duplicated(names(dl))]
# GO annotation system
gsea_res <- dsea_GSEA(drugList=dl, type="GO", ont="MF", exponent=1, nPerm=1000, pAdjustMethod="BH", pvalueCutoff=0.2, minGSSize=5)
gsea_res
result(gsea_res)
# KEGG annotation system
gsea_k_res <- dsea_GSEA(drugList=dl, type="KEGG", exponent=1, nPerm=1000, pAdjustMethod="BH", pvalueCutoff=0.5, minGSSize=5)
gsea_k_res
result(gsea_k_res)
```

# Visulization

## Construct drug-target interaction networks in interesting GO categories

Build drug-target networks in top ranking GO categories
```r
dtnetplot(drugs = dup_hyperG_res@drugs, set = "GO:0032041", ont = "MF", desc="NAD-dependent histone deacetylase activity (H3-K14 specific)")
dtnetplot(drugs = dup_hyperG_res@drugs, set = "GO:0051059", ont = "MF", desc="NF-kappaB binding")
dtnetplot(drugs = hyperG_res@drugs, set = "GO:0001106", ont = "MF", desc="RNA polymerase II transcription corepressor activity")
```

## Construct drug-target interaction networks in interested KEGG pathways or defined other gene/protein sets

Build drug-target networks in top ranking KEGG pathways
```r
dtnetplot(drugs = dup_hyperG_k_res@drugs, set = "hsa05034", desc="Alcoholism")
dtnetplot(drugs = dup_hyperG_k_res@drugs, set = "hsa04330", desc="Notch signaling pathway")
dtnetplot(drugs = dup_hyperG_k_res@drugs, set = "hsa04213", desc="Longevity regulating pathway - multiple species")
```