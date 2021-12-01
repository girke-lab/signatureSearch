# Get instance information of lincs and lincs_expr refdb
library(signatureSearch)

lincs_sig_info <- getTreats("lincs")
save(lincs_sig_info, file="../data/lincs_sig_info.rda")

lincs_expr_inst_info <- getTreats("lincs_expr")
save(lincs_expr_inst_info, file="../data/lincs_expr_inst_info.rda")
