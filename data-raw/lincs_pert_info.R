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

#############################################################
## Add FDA phase and some other information from ChEMBL db ##
#############################################################
#
data("lincs_pert_info")
pcid <- unique(lincs_pert_info$pubchem_cid) 
pcid <- pcid[pcid != "NotFound"] # 8013

# break pcid into batches of size 200 to add compounds to 
# ChemMine Tools (https://chemminetools.ucr.edu/), and
# run Drug-Target Search in ChEMBL db by selecting PubChem as input ID, it takes
# a long time to load and search only 200 compounds, doesn't work well for ~8000 compounds
size <- 200
nbatch <- ceiling(length(pcid)/size)
res_dir <- "~/Desktop/lincs_pcid"
unlink(res_dir, recursive = TRUE)
dir.create(res_dir)
for(i in 1:nbatch){
    sub <- pcid[size*(i-1)+1:min(length(pcid), size*i)]
    writeLines(sub, paste0(res_dir, "/lincs_pcid", i, ".txt"))
}

# Write function to directly query compound annotations from downloaded 
# ChEMBL SQLite db (version 28 when downloading on May 19, 2021) at
# https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/.
# 
library(drugTargetInteractions)
chembldb <- "~/insync/project/ChEMBL_data/chembl_28/chembl_28_sqlite/chembl_28.db"
resultsPath <- "~/insync/project/discover_paper_analysis/results"
config <- genConfig(chemblDbPath=chembldb, resultsPath=resultsPath)
#downloadUniChem(config=config)
#cmpIdMapping(config=config)
#idmap <- readRDS(file.path(config$resultsPath,"cmp_ids.rds"))
data("lincs_pert_info")
pcid <- unique(lincs_pert_info$pubchem_cid); pcid <- pcid[pcid != "NotFound"] # 8013
queryBy <- list(molType="cmp", idType="PubChem_ID", ids=pcid)
qresult <- drugAnnot(queryBy, config=config)
library(magrittr); library(readr)
lincs_pert_info %<>% left_join(qresult, by=c("pubchem_cid"="QueryIDs"))
lincs_pert_info %>% print(n=20, width=Inf)
usethis::use_data(lincs_pert_info, overwrite=TRUE)
write_tsv(lincs_pert_info, "~/insync/project/discover_paper_analysis/data/lincs_pert_info.tsv")

# Function to query known drug annotations in ChEMBL db by providing drug
# "chembl_id", "PubChem_ID" or "DrugBank_ID". 
# The result annotation table contains the following columns:
# 
# "Query ID", "ChEMBL ID", "molregno", "PubChem_ID", "DrugBank_ID", "Preferred Name",
# "Max Phase", "Therapeutic Flag", "Molecule Type",  "First Approval", "Oral", "Parenteral",
# "Topical", "Natural Product", "Inorganic Flag", "USAN Year", "Availability Type", 
# "USAN Stem", "USAN Stem Definition", "Indication Class", "Withdrawn Flag", 
# "Withdrawn Year", "Withdrawn Country", "Withdrawn Reason", "Withdrawn Class",
# "mechanism_of_action", "Action Type", "Direct Interaction", "Molecular Mechanism",
# "Disease Efficacy", "Mechanism Comment", "Selectivity Comment".
# 
drugAnnot <- function(queryBy = list(molType = NULL, idType = NULL, ids = NULL),
                            cmpid_file = file.path(config$resultsPath, "cmp_ids.rds"),
                            config = genConfig()) {
    if (any(names(queryBy) != c("molType", "idType", "ids"))) {
        stop(
            "All three list components in 'queryBy' (named: 'molType',",
            " 'idType' and 'ids') need to be present."
        )
    }
    if (any(vapply(queryBy, length, integer(1)) == 0)) {
        stop(
            "All components in 'queryBy' list need to be populated with ",
            "corresponding character vectors."
        )
    }
    
    ## Load ChEMBL SQLite
    ## ChEMBL SQLite downloaded from here: ftp://ftp.ebi.ac.uk/pub/databases/
    ##    chembl/ChEMBLdb/latest/chembl_24_1_sqlite.tar.gz
    ## after unpacking you get chembl_xx.db
    dbpath <- config$chemblDbPath
    mydb <- dbConnect(SQLite(), dbpath)
    
    ## Input file for following step was generated by cmpIdMapping()
    cmp_ids <- readRDS(cmpid_file)
    rownames(cmp_ids) <- as.character(cmp_ids$molregno)
    
    ## Query by compounds
    if (queryBy$molType == "cmp") {
        ## ID conversions
        if (queryBy$idType == "molregno") {
            cmpvec <- queryBy$ids
        }
        if (queryBy$idType == "chembl_id") {
            cmp_ids <- cmp_ids[!duplicated(cmp_ids$chembl_id), ]
            cmpvec <- as.character(cmp_ids$molregno)
            names(cmpvec) <- as.character(cmp_ids$chembl_id)
            cmpvec <- cmpvec[queryBy$ids]
        }
        if (queryBy$idType == "PubChem_ID") {
            cmp_ids <- cmp_ids[!duplicated(cmp_ids$PubChem_ID), ]
            cmpvec <- as.character(cmp_ids$molregno)
            names(cmpvec) <- as.character(cmp_ids$PubChem_ID)
            cmpvec <- cmpvec[queryBy$ids]
        }
        if (queryBy$idType == "DrugBank_ID") {
            cmp_ids <- cmp_ids[!duplicated(cmp_ids$DrugBank_ID), ]
            cmpvec <- as.character(cmp_ids$molregno)
            names(cmpvec) <- as.character(cmp_ids$DrugBank_ID)
            cmpvec <- cmpvec[queryBy$ids]
        }
        idvec <- paste0("(\"", paste(cmpvec, collapse = "\", \""), "\")")
        
        # "Query ID",  
       # "Direct Interaction", "Molecular Mechanism",
        # "Disease Efficacy", "Mechanism Comment", "Selectivity Comment".
        query <- paste(
            "SELECT a.molregno, a.chembl_id, a.pref_name, a.max_phase,
                a.therapeutic_flag, a.molecule_type, a.first_approval, a.oral,
                a.parenteral, a.topical, a.natural_product, a.inorganic_flag,
                a.usan_year, a.availability_type, a.usan_stem, a.usan_stem_definition,
                a.indication_class, a.withdrawn_flag, a.withdrawn_year, a.withdrawn_country,
                a.withdrawn_reason, a.withdrawn_class,
                b.mechanism_of_action, b.action_type, b.direct_interaction,
                b.molecular_mechanism, b.disease_efficacy, b.mechanism_comment,
                b.selectivity_comment
                      FROM molecule_dictionary AS a
                      LEFT JOIN drug_mechanism AS b ON a.molregno = b.molregno
                      WHERE a.molregno IN", idvec,
                      "GROUP BY a.molregno
                       ORDER BY a.molregno, a.pref_name"
        )
        myquery <- dbSendQuery(mydb, query)
        activityassays <- dbFetch(myquery)
        dbClearResult(myquery)
    }
    resultDF <- data.frame(cmp_ids[as.character(activityassays$molregno), -c(2, 4)], 
                           activityassays[, -c(1, 2)])
    dbDisconnect(mydb)
    
    ## Remove rows with identical values in all fields
    resultDF <- resultDF[!duplicated(apply(resultDF, 1, paste, collapse = "_")), ]
    ## Add query column and sort rows in result table according to query
    resultDF <- data.frame(
        QueryIDs = names(cmpvec),
        resultDF[match(names(cmpvec), resultDF[[queryBy$idType]]), ]
    )
    rownames(resultDF) <- NULL
    return(tibble::tibble(resultDF))
}
