# Code used to generate lincs_pert_info.rds
# 
###############################
## Fix Broad LINCS pert info ##
###############################
#
### Fix PCID -666 terms
# Some PubChem CID annotation of Broad LINCS compounds in the 
# `GSE92742_Broad_LINCS_pert_info.txt` file frequently show `-666` entries.
# However, most compounds with -666 are present in PubChem. The following code 
# fixes this issue by querying PubChem webpage by compound names to get the 
# corresponding PCID.
source("getPertInfo.R")
lincs_pert_info <- read_tsv("data/GSE92742_Broad_LINCS_pert_info.txt")
pert2 <- lincs_pert_info %>% distinct(pert_iname, .keep_all=TRUE) %>% 
    filter(pert_type == "trt_cp") # 19,811 x 8
pert2[is.na(pert2$pubchem_cid), "pubchem_cid"] <- -666
sum(pert2$pubchem_cid == "-666", na.rm = TRUE) # 916
# Replace -666 if there are actually PCID
pert_nopcid <- pert2$pert_iname[pert2$pubchem_cid=="-666"]
pcid_new <- vapply(pert_nopcid, function(i) name2pubchemCID(i), character(1)) # takes long time
sum(pcid_new=="NotFound") # 84
pert2$pubchem_cid[pert2$pubchem_cid=="-666"] <- pcid_new
write_tsv(pert2, "data/broad_lincs_pert_info.tsv")

# 
##############################################################
## Subset Broad LINCS pert info to compounds in LINCS refdb ##
##############################################################
#
# Since most Broad LINCS pert info are BRD drugs, the non-BRD drugs are almost
# all present in the LINCS refdb, so subset the broad pert info to reduce size
# 
library(readr)
broad <- read_tsv("../discover_paper_analysis/data/broad_lincs_pert_info.tsv")
data("lincs_sig_info")
lincs_pert_info <- broad[broad$pert_iname %in% lincs_sig_info$pert,] # 8140 x 8
usethis::use_data(lincs_pert_info, overwrite=TRUE)
