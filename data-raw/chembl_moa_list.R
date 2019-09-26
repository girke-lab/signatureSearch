## code to prepare `chembl_moa_list` dataset goes here

## MOA to gene entrez id mappings from ChEMBL database

# The MOA categories from ChEMBL, can be obtained from the drug_mechanism table 
# in ChEMBL SQLite database with
library(RSQLite)
mydb <- dbConnect(SQLite(), "~/insync/project/ChEMBL_data/chembl_25/chembl_25_sqlite/chembl_25.db")
# chembl_25.db is downloaded at ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_25_sqlite.tar.gz
drug_mech <- dbGetQuery(mydb, 'SELECT * FROM drug_mechanism')
# in this table (here drug_mech) the mechanism_of_action column contains 1287 unique entries:
length(unique(drug_mech$mechanism_of_action))
# The corresponding target proteins mapping to the mechanism_of_action categories 
# can be obtained via the target_dictionary table with the tid (target id) column as common key
target_dict <- dbGetQuery(mydb, 'SELECT * FROM target_dictionary')
# Just as a reminder, the tables and column definitions from ChEMBL are defined in 
# the *_schema_documentation.html file here: ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/

# The target ids from ChEMBL are given as ChEMBL protein ids. Those can be translated 
# into any other protein id type (e.g. UniProt) via UniChem.
moa_target <- unique(inner_join(drug_mech[,c("mechanism_of_action", "tid", "action_type")],
                                target_dict[,c("tid","target_type", "pref_name","tax_id","organism","chembl_id")], by="tid"))
# Save chembl_id in the 'moa_target' to a csv file used for pyhton ChEMBL webresource client
# to map to Uniprot id.
write_csv(data.frame(chembl_id=moa_target[,'chembl_id']), "data/chembl_moa_target_chembl_id.csv")
# Read the target chembl id to uniprot id mappings generated from "target_chembl2uniprotID.ipynb"
chem2uni <- read_csv("data/chembl_moa_target_chem2uni.csv", col_names=FALSE)
colnames(chem2uni) = c("chembl_id", "uniprot_id")
moa_target %<>% inner_join(chem2uni, by="chembl_id")
# Get moa to uniprot id mappings
moa2uni <- unique(na.omit(moa_target[,c("mechanism_of_action", "uniprot_id")]))
## Combine uniprot ids of duplicated moa
moa2uni_aggr <- aggregate(.~mechanism_of_action, data=moa2uni, FUN=paste0, collapse="; ")
moa_uni_list <- lapply(moa2uni_aggr$uniprot_id, function(x) unique(unlist(strsplit(x, '; '))))
names(moa_uni_list) <- moa2uni_aggr$mechanism_of_action
moa_uni_df <- data.frame("moa"=rep(names(moa_uni_list), times=sapply(moa_uni_list, length)),
                         "UNIPROT"=unlist(moa_uni_list), stringsAsFactors=FALSE)
## Convert uniprot id to entrez id in homo sapiens
uni_entr_df = select(org.Hs.eg.db, keys=unique(moa_uni_df$UNIPROT), 
                     columns=c("ENTREZID","UNIPROT"), keytype="UNIPROT")
moa_uni_entr <- inner_join(moa_uni_df, uni_entr_df, by="UNIPROT")
moa_entr_df <- unique(na.omit(moa_uni_entr[,c("moa","ENTREZID")]))
write_csv(moa_entr_df, "data/ChEMBL_moa2entrezid_df.csv")
moa_entr_list <- split(moa_entr_df$ENTREZID, moa_entr_df$moa)
saveRDS(moa_entr_list, "data/ChEMBL_moa2entrezid_list.rds")
chembl_moa_list <- readRDS("data/ChEMBL_moa2entrezid_list.rds")

usethis::use_data("chembl_moa_list")
