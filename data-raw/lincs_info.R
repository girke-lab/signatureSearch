# Get instance information of lincs and lincs_expr refdb
library(signatureSearch)s

lincs_sig_info <- get_treat_info("lincs")
save(lincs_sig_info, file="../data/lincs_sig_info.rda")

lincs_expr_inst_info <- get_treat_info("lincs_expr")
save(lincs_expr_inst_info, file="../data/lincs_expr_inst_info.rda")
