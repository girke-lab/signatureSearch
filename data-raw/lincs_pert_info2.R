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
write_tsv(cmpinfo2, "~/insync/project/discover_paper_analysis/data/lincs_pert_info2.tsv")
cmpinfo2 <- read_tsv("~/insync/project/discover_paper_analysis/data/lincs_pert_info2.tsv")

##################################
## Add pubchem_cid and many other annotations from 2017 lincs_pert_info 
#################################
data("lincs_pert_info")
cmpinfo2 %<>% left_join(select(
    lincs_pert_info, -c("pert_id", "pert_type", "MOA", "inchi_key_prefix", "inchi_key", "canonical_smiles")), 
    by="pert_iname")
write_tsv(cmpinfo2, "~/insync/project/discover_paper_analysis/data/lincs_pert_info2.tsv")
lincs_pert_info2 <- cmpinfo2
usethis::use_data(lincs_pert_info2, overwrite=TRUE)
