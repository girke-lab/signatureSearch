## code to prepare `clue_moa_list` dataset goes here

## GO to QUERY tool in CLUE website (https://clue.io/query), use example query 
## signature or custom one (e.g. Peters) to query the Touchstone database (Version: 1.0.1.1).
## Go to Query History to download the result txt file to the local under data 
## direcotry as "conn.txt"
library(readr)
source("getPertInfo.R")
conn_res <- read_tsv("data/conn.txt")
lincs_drugs <- unique(as.character(conn_res[,"Name"]))
#lincs_drugs <- readLines("~/insync/project/longevityTools_eDRUG/data/lincs_drugs.txt") # 2512

## Access to the LINCS/CLUE api and clue.io requires a user key, which can be obtained 
## from the LINCS by registering at https://clue.io (free for non-commercial use). 
lincs_pert_info <- getPertInfo(lincs_drugs, user_key="yourUserKey") # Take some time
saveRDS(lincs_pert_info, "data/lincs_pert_info_api.rds")
#lincs_pert_info <- readRDS("~/insync/project/lincs_gse92742_dataset_analysis/data/lincs_pert_info_api.rds")

## add target number column to lincs_pert_info
rep_num <- sapply(lincs_pert_info$moa, length)
pert_moa <- na.omit(unique(data.frame(pert_iname=rep(lincs_pert_info$pert_iname, rep_num),  
                                      moa=unlist(lincs_pert_info$moa), stringsAsFactors = FALSE)))
moa_list <- split(pert_moa$pert_iname, pert_moa$moa)
idx <- sapply(moa_list, length)>1
moa_list <- moa_list[idx] # 346
clue_moa_list <- moa_list[names(moa_list)!="-666"] #345 
usethis::use_data("clue_moa_list")
