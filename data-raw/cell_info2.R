## code to prepare `cell_info2` dataset goes here
setwd("~/insync/project/discover_paper_analysis/")
# download.file("https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/cellinfo_beta.txt", 
#              "downloads/cellinfo_beta.txt")
library(data.table)
cellInfo2020 <- fread("downloads/cellinfo_beta.txt")
head(cellInfo2020, 20)

#### Load source data for cell_info file from 2017 version ###
# download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&format=file&file=GSE92742%5FBroad%5FLINCS%5Fcell%5Finfo%2Etxt%2Egz", 
#              "downloads/GSE92742_Broad_LINCS_cell_info.txt.gz")
cellInfo2017 <- fread("downloads/GSE92742_Broad_LINCS_cell_info.txt.gz")
head(cellInfo2017)

#### create new "cell_info" table ####
setnames(cellInfo2020, "cell_iname", "cell_id")
cell_info2020 <- merge(cellInfo2020, cellInfo2017[,c("cell_id", "primary_site")], 
                       by = "cell_id", all.x = TRUE)
dim(cell_info2020)
head(cell_info2020)
cell_info2020[160, "subtype"] <- "medullary cystic kidney disease (MCKD) type 1"
library(readr)
write_tsv(cell_info2020, "data/cell_info2.tsv")

cell_info2 <- read_tsv("~/insync/project/discover_paper_analysis/data/cell_info2.tsv")
usethis::use_data(cell_info2, overwrite=TRUE)
