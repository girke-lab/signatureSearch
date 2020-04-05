#' Run the entire GESS/FEA workflow
#' 
#' This function runs the entire GESS/FEA workflow when providing the query
#' drug and cell type, as well as selecting the reference database (e.g.
#' cmap or lincs), defining the specific GESS and FEA methods. In this case,
#' the query GES is drawn from the reference database. The GESS/FEA 
#' results will be a list object in R session. It will also writes the GESS/FEA 
#' result table to user's current working directory by making a working 
#' environment (creating a directory named by use case). The result tables are 
#' stored in the \code{results} directory, the working environment also contains
#' a template Rmd vignette, users could make some modification according to 
#' their use case and render it to generate an HTML report.
#' 
#' @param drug character(1) representing query drug name (e.g. vorinostat). 
#' This query drug should be included in the \code{refdb}
#' @param cell character(1) indicating the cell type that the query drug
#' treated in
#' @param refdb character(1), one of "lincs", "lincs_expr", "cmap", "cmap_expr",
#' or path to the HDF5 file built from \code{\link{build_custom_db}} function
#' @param gess_method character(1), one of "LINCS", "CORsub", "CORall", 
#' "Fisher", "CMAP", "gCMAP"
#' @param fea_method character(1), one of "dup_hyperG", "mGSEA", "mabs", 
#' "hyperG", "GSEA"
#' @param N_gess_drugs number of unique drugs in GESS result used as input of
#' FEA 
#' @param Nup integer(1). Number of most up-regulated genes to be subsetted 
#' for GESS query when \code{gess_method} is CMAP or LINCS
#' @param Ndown integer(1). Number of most down-regulated genes to be subsetted 
#' for GESS query when \code{gess_method} is CMAP or LINCS
#' @param higher numeric(1), it is defined when gess_method argument is 'gCMAP'
#' or 'Fisher' representing the 'upper' threshold of subsetting genes with a 
#' score larger than 'higher'
#' @param lower numeric(1), it is defined when gess_method argument is 'gCMAP'
#' or 'Fisher' representing the 'lower' threshold of subsetting genes
#' @return 
#' @examples 
#' drug="vorinostat"; cell="SKB"
#' refdb <- system.file("extdata", "sample_db.h5", package="signatureSearch")
runWF <- function(drug, cell, refdb, gess_method, fea_method, N_gess_drugs,
                  Nup=150, Ndown=150, higher=1, lower=-1){
    if(gess_method == "CMAP"){
        query <- getDEGSig(drug, cell, refdb=refdb)
        qsig_cmap <- qSig(query = list(upset=query[[1]], downset=query[[2]]), 
                          gess_method="CMAP", refdb=refdb)
        gess_res <- gess_cmap(qSig=qsig_cmap, chunk_size=5000)
    }  
    if(gess_method == "LINCS"){
        query <- getDEGSig(drug, cell, refdb=refdb)
        qsig_lincs <- qSig(query = list(upset=query[[1]], downset=query[[2]]), 
                           gess_method="LINCS", refdb=refdb)
        gess_res <- gess_lincs(qsig_lincs, sortby="NCS", tau=TRUE)
    }
    if(gess_method == "gCMAP"){
        query <- getSig(drug, cell, refdb=refdb)
        qsig_gcmap <- qSig(query=query, gess_method="gCMAP", refdb=refdb)
        gess_res <- gess_gcmap(qsig_gcmap, higher=higher, lower=lower)
    }
    if(gess_method == "Fisher"){
        query <- getSig(drug, cell, refdb=refdb)
        qsig_fisher <- qSig(query=query, gess_method="Fisher", refdb=refdb)
        gess_res <- gess_fisher(qSig=qsig_fisher, higher=higher, lower=lower)
    }
    
    qsig_sp <- qSig(query = query_mat, gess_method = "Cor", refdb = db_path)
    sp <- gess_cor(qSig=qsig_sp, method="spearman")
    result(sp)
    
    query_mat_sub <- as.matrix(query_mat[c(upset, downset),])
    qsig_spsub <- qSig(query = query_mat_sub, gess_method = "Cor", refdb = db_path)
    spsub <- gess_cor(qSig=qsig_spsub, method="spearman")
    result(spsub)
    
    drugs <- unique(result(lincs)$pert[1:10])
    dup_hyperG_res <- tsea_dup_hyperG(drugs = drugs, universe = "Default", 
                                      type = "GO", ont="MF", pvalueCutoff=0.05, 
                                      pAdjustMethod="BH", qvalueCutoff = 0.1, 
                                      minGSSize = 10, maxGSSize = 500)
    dup_hyperG_res
    
    dup_hyperG_k_res <- tsea_dup_hyperG(drugs = drugs, universe = "Default", 
                                        type = "KEGG", pvalueCutoff=0.5, 
                                        pAdjustMethod="BH", qvalueCutoff = 0.5, 
                                        minGSSize = 10, maxGSSize = 500)
    result(dup_hyperG_k_res)
    
    mgsea_res <- tsea_mGSEA(drugs=drugs, type="GO", ont="MF", exponent=1, 
                            nPerm=1000, pvalueCutoff=1, minGSSize=5)
    result(mgsea_res)
    
    mgsea_k_res <- tsea_mGSEA(drugs=drugs, type="KEGG", exponent=1, 
                              nPerm=1000, pvalueCutoff=1, minGSSize=2)
    result(mgsea_k_res)
    
    mabs_res <- tsea_mabs(drugs=drugs, type="GO", ont="MF", nPerm=1000, 
                          pvalueCutoff=0.05, minGSSize=5)
    result(mabs_res)
    
    mabs_k_res <- tsea_mabs(drugs=drugs, type="KEGG", nPerm=1000, 
                            pvalueCutoff=0.2, minGSSize=5)
    result(mabs_k_res)
    
    drugs <- unique(result(lincs)$pert[1:10])
    hyperG_res <- dsea_hyperG(drugs=drugs, type="GO", ont="MF")
    result(hyperG_res)
    
    hyperG_k_res <- dsea_hyperG(drugs = drugs, type = "KEGG", 
                                pvalueCutoff = 1, qvalueCutoff = 1, 
                                minGSSize = 10, maxGSSize = 2000)
    result(hyperG_k_res)
    
    dl <- abs(result(lincs)$NCS); names(dl) <- result(lincs)$pert
    dl <- dl[dl>0]
    dl <- dl[!duplicated(names(dl))]
    gsea_res <- dsea_GSEA(drugList=dl, type="GO", ont="MF", exponent=1, nPerm=1000,
                          pvalueCutoff=0.2, minGSSize=5)
    result(gsea_res)
    
    gsea_k_res <- dsea_GSEA(drugList=dl, type="KEGG", exponent=1, nPerm=1000, 
                            pvalueCutoff=1, minGSSize=5)
    result(gsea_k_res)
    
} 