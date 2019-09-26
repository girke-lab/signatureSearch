## code to prepare `cell_info` dataset goes here

## Download and unzip the GSE92742_Broad_LINCS_cell_info.txt.gz file
## from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
cell_info <- read_tsv("~/insync/project/lincs_gse92742_dataset_analysis/data/GSE92742_Broad_LINCS_cell_info.txt")
cell_name <- c("NEU", "A549", "PC3", "VCAP", "MCF7", "ASC", "SKB", "PHH", "NPC",
               "HT29", "FIBRNPC", "HCC515", "A375", "HA1E", "CD34", "HL60", "U937", "HEK293T",
               "HEPG2", "HUH7", "NOMO1", "THP1", "BT20", "HS578T", "MCF10A", "MDAMB231", "SKBR3",   
               "NKDBA", "JURKAT", "U266")
# 30 cell types tested in `lincs` database from the `signatureSearchData` package  
cell_info %<>% dplyr::filter(base_cell_id %in% cell_name & 
                !duplicated(paste(cell_info$base_cell_id, cell_info$sample_type, sep="_"))) %>%
    dplyr::select(c("cell_id","sample_type","primary_site","subtype")) %>% 
    dplyr::rename("cell_type"="sample_type")
# replace "primary" as "normal" in sample_type column
cell_info %<>% mutate(cell_type=gsub("primary","normal",cell_type))
#write_tsv(cell_info, "~/insync/project/GESS_and_FEA/data/cell_info.tsv")
usethis::use_data("cell_info")
