# Get instance information of lincs_expr refdb
library(signatureSearch); library(HDF5Array); library(ExperimentHub)
eh <- ExperimentHub()
lincs_expr_path <- eh[["EH3227"]]
inst <- as.character(HDF5Array(lincs_expr_path, name="colnames")) # 38824
inst_df <- as.data.frame(t(vapply(seq_along(inst), function(i)
    unlist(strsplit(as.character(inst[i]), "__")), 
    FUN.VALUE=character(3))), stringsAsFactors=FALSE)
colnames(inst_df) <- c("pert", "cell", "pert_type")
lincs_expr_inst_info <- as_tibble(inst_df)
save(lincs_expr_inst_info, file="../data/lincs_expr_inst_info.rda")
