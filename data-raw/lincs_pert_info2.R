# Code used to generate lincs_pert_info2.rda

# lincs_pert_info2 contains compound annotations for the new LINCS 2020 beta 
# database at https://clue.io/releases/data-dashboard

########################################################
## Retrieve LINCS beta 2020 compound annotation table ##
########################################################

# download.file("https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/compoundinfo_beta.txt",
#             "~/insync/project/discover_paper_analysis/data/compoundinfo_beta.txt")
cmpinfo <- read_tsv("~/insync/project/discover_paper_analysis/data/compoundinfo_beta.txt") # 39321 x 7

# Issue: Some pert_id have different cmap_name.
id2name <- split(cmpinfo$cmap_name, cmpinfo$pert_id)
id2name <- lapply(id2name, unique)
index <- which(sapply(id2name, length) > 1) # 18
id2name[index]
# Solution: get the first cmap_name by unique
length(unique(cmpinfo$pert_id)) # 34419

## Collapse multiple targets and MOAs for pert_id by '; '
id2tar <- split(cmpinfo$target, cmpinfo$pert_id)
id2tar_vec <- sapply(id2tar, function(x) paste(unique(x), collapse="; "))
id2moa <- split(cmpinfo$moa, cmpinfo$pert_id)
id2moa_vec <- sapply(id2moa, function(x) paste(unique(x), collapse="; "))
cmpinfo2 <- distinct(cmpinfo[,c("pert_id", "cmap_name", "canonical_smiles", 
                                "inchi_key", "compound_aliases")], 
                     pert_id, .keep_all=TRUE)
cmpinfo2 %<>% mutate(target=id2tar_vec[cmpinfo2$pert_id], MOA=id2moa_vec[cmpinfo2$pert_id]) %>%
    dplyr::rename(pert_iname=cmap_name) # 34419 x 7

##################################
## Add pubchem_cid and many other annotations from 2017 lincs_pert_info 
#################################
data("lincs_pert_info")
cmpinfo2 %<>% left_join(select(
    lincs_pert_info, -c("pert_id", "pert_type", "MOA", "inchi_key_prefix", "inchi_key", "canonical_smiles")), 
    by="pert_iname")

##########################################################
## Including MOAclue, mergeMOA, t_gn_sym, mergeTargets, Target_pathway columns
##########################################################
setwd("~/insync/project/discover_paper_analysis")

## mergeMOA: MOA merged from LINCS2 compoundinfo_beta.txt, CLUE moa list and 
## mechanism_of_action from ChEMBL.
## mergeTargets: merged targets from t_gn_sym and Targets column from compoundinfo_beta.txt

#### merge pert_iname and compound_aliases column together ####
name <- cmpinfo2$pert_iname; alias <- cmpinfo2$compound_aliases; NameAliasMer <- NULL
for(i in 1:length(name)){
    if(!grepl("BRD-", name[i])){
        NameAliasMer[i] <- name[i]
    } else {
        if(is.na(alias[i])){
            NameAliasMer[i] <- name[i]
        } else {
            NameAliasMer[i] <- alias[i]
        }
    }
}
library(dplyr); library(readr); library(magrittr); library(data.table)
cmpinfo2 %<>% mutate(NameAliasMer=NameAliasMer)
finalDT <- cmpinfo2

#### Add MOAclue ####
addMOA <- function(df, drug_col, drug2moa_path, moa_col="MOA"){
    drug2moa <- readRDS(drug2moa_path)
    drug2moa_clp <- lapply(drug2moa, paste, collapse="; ")
    d2m_vec <- unlist(drug2moa_clp); names(d2m_vec) <- names(drug2moa_clp)
    df[[moa_col]] <- d2m_vec[df[[drug_col]]]
    return(df)}
finalDT <- addMOA(finalDT, drug_col="NameAliasMer", moa_col="MOAclue",
                  drug2moa_path="data/clue_ts_drug2moa_list_1.2.rds")

#### Merge MOA, MOAclue, mechanism_of_action columns together ####
MOA1 <- finalDT$MOA; MOA2 <- finalDT$MOAclue; MOA3 <- finalDT$mechanism_of_action;
mergeMOA <- NULL
for(i in 1:length(MOA1)){
    merge <- paste(na.omit(c(MOA1[i], MOA2[i], MOA3[i])), collapse="; ")
    moa <- unique(unlist(strsplit(merge, split="; ")))
    mergeMOA[i] <- paste(moa, collapse="; ")
}

finalDT %<>% mutate(mergeMOA = mergeMOA)

#### add t_gn_sym, mergeTargets, Ntar information ####
tars <- get_targets(finalDT$NameAliasMer)
finalDT %<>% left_join(tars, by=c("NameAliasMer"="drug_name"))

target1 <- finalDT$target; target2 <- finalDT$t_gn_sym
mergeTargets <- NULL
for(i in 1:length(target1)){
    merge <- paste(na.omit(c(target1[i], target2[i])), collapse="; ")
    tar <- unique(unlist(strsplit(merge, split="; ")))
    mergeTargets[i] <- paste(tar, collapse="; ")
}
finalDT %<>% mutate(mergeTargets=mergeTargets)
finalDT$Ntar <- sapply(finalDT$mergeTargets, function(tar) {
    if(is.na(tar)) return(0); length(unlist(strsplit(tar, split="; ")))})


#### Add LAD information ####
#############################

######################## Get overlapped drugs in LINCS2 and DrugAge4 (run once)
da_annot <- read_csv("downloads/drugage_clean.csv")
writeLines(da_annot$compound_name, "data/drugs_drugage4.txt")
writeLines(finalDT$NameAliasMer, "data/nameAliasMer_lincs2.txt")
# Get drugs common in the newest DrugAge database (Build 4 (20/11/2021)) and 
# LINCS2 database by PubChem CID matching. Specifically, using the 
# 'data/drugs_drugage4.txt' and 'data/nameAliasMer_lincs2.txt' file as input 
# Synonyms to the PubChem Identifier
# Exchange Service at https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi.
# Operator Type is selected as Same CID, Output IDs is selected as CIDs. Output
# Method is selected as Two column file. Compression is selected as No compression.
# Then Submit Job and save the result to the 'downloads/drugs_da4_name2pcid.txt' and
# 'downloads/nameAliasMer_lincs2_name2pcid.txt' file
da_sym2cid <- read_tsv("downloads/drugs_da4_name2pcid.txt", col_names=FALSE)
lincs_sym2cid <- read_tsv("downloads/nameAliasMer_lincs2_name2pcid.txt", col_names=FALSE)

# common drugs in LINCS and DrugAge database
cmcid <- na.omit(unique(intersect(da_sym2cid$X2, lincs_sym2cid$X2)))
length(cmcid) # 375
drug_cm <- unique(lincs_sym2cid$X1[lincs_sym2cid$X2 %in% cmcid]) # 253
writeLines(drug_cm, "data/drugs_lincs2_and_drugage4.txt")
#############################################

lincs_da_drugs <- readLines("data/drugs_lincs2_and_drugage4.txt")
#### annotate by pubchem_cid ####
finalDT$isLAD <- finalDT$NameAliasMer %in% lincs_da_drugs

#### annotate with Reactome pathways target is in ####
# gmtpath is at "/rhome/bgongol/bigdata/longevitySignatures_RNAseq/data/reactome/ReactomePathways.gmt"
# on cluster, copied under downloads
paths <- clusterProfiler::read.gmt("downloads/ReactomePathways.gmt")
length(unique(paths$term)) # 2504
targs <- finalDT$mergeTargets
TargetPathwayMatch <- function(targets, pathways){
    PathwayID <- NULL
    for(i in 1:length(targets)){
        spl <- unlist(strsplit(targets[i], split="; "))
        if(length(spl)==0){
            PathwayID[i] <- ""
        } else {
            IDd <- unique(pathways[(pathways$gene %in% spl),]$term)
            CombPath <- paste(IDd, collapse = "; ")
            PathwayID[i] <- CombPath}}
    return(PathwayID)
}
PathwayID <- TargetPathwayMatch(targets=targs, pathways=paths)
finalDT$Target_pathway <- PathwayID

lincs_pert_info2 <- finalDT # 34419 x 48
write_tsv(lincs_pert_info2, "~/insync/project/discover_paper_analysis/data/lincs_pert_info2.tsv")
lincs_pert_info2 <- read_tsv("~/insync/project/discover_paper_analysis/data/lincs_pert_info2.tsv")
usethis::use_data(lincs_pert_info2, overwrite=TRUE)


