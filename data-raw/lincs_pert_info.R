# Code used to generate lincs_pert_info.rda
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
getwd() # start from the signatureSearch directory
source("data-raw/getPertInfo.R")
# The GSE92742_Broad_LINCS_pert_info.txt file was downloaded and unzipped 
# from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&format=file&file=GSE92742%5FBroad%5FLINCS%5Fpert%5Finfo%2Etxt%2Egz 
# to the data folder.
library(readr); library(dplyr); library(magrittr)
lincs_pert_info <- read_tsv("~/insync/project/discover_paper_analysis/data/GSE92742_Broad_LINCS_pert_info.txt")
pert2 <- lincs_pert_info %>% distinct(pert_iname, .keep_all=TRUE) %>% 
    filter(pert_type == "trt_cp") # 19,811 x 8
pert2[is.na(pert2$pubchem_cid), "pubchem_cid"] <- -666
sum(pert2$pubchem_cid == "-666", na.rm = TRUE) # 916
# Replace -666 if there are actually PCID
pert_nopcid <- pert2$pert_iname[pert2$pubchem_cid=="-666"]
pcid_new <- vapply(pert_nopcid, function(i) name2pubchemCID(i), character(1)) # takes about 10 minutes
sum(pcid_new=="NotFound") # 84
pert2$pubchem_cid[pert2$pubchem_cid=="-666"] <- pcid_new
write_tsv(pert2, "~/insync/project/discover_paper_analysis/data/broad_lincs_pert_info_fixed.tsv")

##############################################################
## Subset Broad LINCS pert info to compounds in LINCS refdb ##
##############################################################

# Since most Broad LINCS pert info are BRD drugs, the non-BRD drugs are almost
# all present in the LINCS refdb, so subset the broad pert info to reduce size
broad <- read_tsv("~/insync/project/discover_paper_analysis/data/broad_lincs_pert_info_fixed.tsv")
data("lincs_sig_info")
lincs_pert_info <- broad[broad$pert_iname %in% lincs_sig_info$pert,] # 8140 x 8

########################
## Add MOA annotation ##
########################
# Add drug MOA from CLUE Touchstone compound annotations
data("clue_moa_list")
lincs_pert_info <- addMOA(lincs_pert_info, "pert_iname", moa_list=clue_moa_list)

#############################################################
## Add FDA phase and some other information from ChEMBL db ##
#############################################################
pcid <- unique(lincs_pert_info$pubchem_cid) 
pcid <- pcid[pcid != "NotFound"] # 8013

# Use PubChem CID to query compound annotations from downloaded 
# ChEMBL SQLite db (version 28 when downloading on May 19, 2021) at
# https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/.
# 
library(drugTargetInteractions)
chembldb <- "~/insync/project/ChEMBL_data/chembl_28/chembl_28_sqlite/chembl_28.db"
resultsPath <- "~/insync/project/discover_paper_analysis/results"
config <- genConfig(chemblDbPath=chembldb, resultsPath=resultsPath)
#downloadUniChem(config=config)
cmpIdMapping(config=config)
#idmap <- readRDS(file.path(config$resultsPath,"cmp_ids.rds"))
pcid <- unique(lincs_pert_info$pubchem_cid); pcid <- pcid[pcid != "NotFound"] # 8013
queryBy <- list(molType="cmp", idType="PubChem_ID", ids=pcid)
qresult <- drugAnnot(queryBy, config=config)

library(magrittr); library(readr)
lincs_pert_info %<>% left_join(qresult, by=c("pubchem_cid"="QueryIDs"))
lincs_pert_info %>% print(n=20, width=Inf)
usethis::use_data(lincs_pert_info, overwrite=TRUE)
write_tsv(lincs_pert_info, "~/insync/project/discover_paper_analysis/data/lincs_pert_info.tsv")


