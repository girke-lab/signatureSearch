## code to prepare `cell_info2` dataset goes here
# cell_info2 <- fread("/home/bgongol/Downloads/Updatefiles/cell_info2.xls")
#library(usethis)
#use_data_raw()

library(data.table)
download.file("https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/cellinfo_beta.txt",
              "/home/bgongol/Downloads/Updatefiles/cellinfo_beta.txt")
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz",
              "/home/bgongol/Downloads/Updatefiles/GSE92742_Broad_LINCS_cell_info.txt.gz")
system("gunzip /home/bgongol/Downloads/Updatefiles/GSE92742_Broad_LINCS_cell_info.txt.gz")

LINCSCellInfoGen <- function(Fpath2 = cellInfo2Path, Fpath2017 = cellInfo2017Path){
  cellInfo2 <- fread(Fpath2)
  cell_info2 <- data.frame(cell_id = cellInfo2$cell_iname,
                           cell_type = cellInfo2$cell_type,
                           # primary_site = cellInfo2$
                           subtype = cellInfo2$subtype,
                           donor_sex = cellInfo2$donor_sex)
  cellInfo2017 <- fread(Fpath2017)
  cell_info2 <- merge(cell_info2, cellInfo2017[,c("cell_id", "primary_site")], by = "cell_id", all.x = TRUE)
  return(cell_info2)}

cell_info2 <- LINCSCellInfoGen(Fpath2 = "/home/bgongol/Downloads/Updatefiles/cellinfo_beta.txt",
                 Fpath2017 = "/home/bgongol/Downloads/Updatefiles/GSE92742_Broad_LINCS_cell_info.txt")
cell_info2

usethis::use_data(cell_info2, overwrite = TRUE)
