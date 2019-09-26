## code to prepare `univ_go` dataset

library(org.Hs.eg.db)
go <- clusterProfiler:::get_GO_data(OrgDb = org.Hs.eg.db, ont="ALL", keytype = "SYMBOL")
names(go)
goterms <- get("PATHID2EXTID", go)
genes <- unique(unlist(goterms)) # 19456
# The length of `genes` may be different from 19456 due to update of the `org.Hs.eg.db` package
univ_go <- genes
#writeLines(genes, "~/insync/project/GESS_and_FEA/data/univ_genes_in_GO.txt")
use_data(univ_go, internal = TRUE)
