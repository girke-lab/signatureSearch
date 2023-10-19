#' Run the Entire GESS/FEA Workflow
#' 
#' @description This function runs the entire GESS/FEA workflow when providing 
#' the query drug and cell type, as well as selecting the reference database 
#' (e.g. 'cmap' or 'lincs'), defining the specific GESS and FEA methods. 
#' In this case, the query GES is drawn from the reference database. 
#' The N (defined by the `N_gess_drugs` argument) top ranking hits in the GESS 
#' tables were then used for FEA where three different annotation systems were 
#' used: GO Molecular Function (GO MF), GO Biological Process (GO BP) and 
#' KEGG pathways. 
#' 
#' The GESS/FEA results will be stored in a list object in R session. 
#' A working environment named by the use case will be created under users 
#' current working directory or under other directory defined by users.
#' This environment contains a \code{results} folder where the GESS/FEA 
#' result tables were written to. The working environment also contains
#' a template Rmd vignette as well as a rended HTML report, users could make
#' modifications on the Rmd vignette as they need and re-render it to generate 
#' their HTML report.
#' 
#' @param drug character(1) representing query drug name (e.g. vorinostat). 
#' This query drug should be included in the \code{refdb}
#' @param cell character(1) indicating the cell type that the query drug
#' treated in. Details about cell type options in LINCS database can be found 
#' in the \code{cell_info} table after load the `signatureSearch` package
#' and running `data("cell_info")`
#' @param refdb character(1), one of "lincs", "lincs_expr", "cmap", "cmap_expr",
#' or path to the HDF5 file built from \code{\link{build_custom_db}} function
#' @param gess_method character(1), one of "LINCS", "CORsub", "CORall", 
#' "Fisher", "CMAP", "gCMAP". When \code{gess_method} is "CORsub" or "CORall",
#' only "lincs_expr" or "cmap_expr" databases are supported.
#' @param fea_method character(1), one of "dup_hyperG", "mGSEA", "mabs", 
#' "hyperG", "GSEA"
#' @param N_gess_drugs number of unique drugs in GESS result used as input of
#' FEA 
#' @param env_dir character(1), directory under which the result environment 
#' located. The default is users current working directory in R session, can
#' be checked via \code{getwd()} command in R  
#' @param tau TRUE or FALSE indicating whether to compute Tau scores if 
#' \code{gess_method} is set as 'LINCS'
#' @param Nup integer(1). Number of most up-regulated genes to be subsetted 
#' for GESS query when \code{gess_method} is CMAP, LINCS or CORsub
#' @param Ndown integer(1). Number of most down-regulated genes to be subsetted 
#' for GESS query when \code{gess_method} is CMAP, LINCS or CORsub
#' @param higher numeric(1), it is defined when gess_method argument is 'gCMAP'
#' or 'Fisher' representing the 'upper' threshold of subsetting genes with a 
#' score larger than 'higher'
#' @param lower numeric(1), it is defined when gess_method argument is 'gCMAP'
#' or 'Fisher' representing the 'lower' threshold of subsetting genes
#' @param method One of 'spearman' (default), 'kendall', or 'pearson', 
#' indicating which correlation coefficient to use
#' @param pvalueCutoff double, p-value cutoff for FEA result
#' @param qvalueCutoff double, qvalue cutoff for FEA result
#' @param minGSSize integer, minimum size of each gene set in annotation system
#' @param maxGSSize integer, maximum size of each gene set in annotation system
#' @param runFEA Logical value indicating if FEA analysis is performed.
#' @param GenerateReport Logical value indicating if a report is generated.
#' @return list object containing GESS/FEA result tables
#' @importFrom readr write_tsv
#' @export
#' @examples 
#' drug <- "vorinostat"; cell <- "SKB"
#' refdb <- system.file("extdata", "sample_db.h5", package="signatureSearch")
#' env_dir <- tempdir()
#' wf_list <- runWF(drug, cell, refdb, gess_method="LINCS", 
#'     fea_method="dup_hyperG", N_gess_drugs=10, env_dir=env_dir, tau=FALSE,
#'     runFEA=FALSE, GenerateReport= FALSE)
#'
runWF <- function(drug, cell, refdb, gess_method, fea_method, N_gess_drugs=100,
                  env_dir=".", tau=TRUE, Nup=150, Ndown=150, higher=1, lower=-1, 
                  method="spearman", pvalueCutoff=1, qvalueCutoff=1, minGSSize=5, 
                  maxGSSize=500, runFEA=TRUE, GenerateReport= TRUE){
  if(gess_method == "CMAP"){
    query <- getDEGSig(drug, cell, Nup=Nup, Ndown=Ndown, refdb=refdb)
    qsig_cmap <- qSig(query = list(upset=query[[1]], downset=query[[2]]), 
                      gess_method="CMAP", refdb=refdb)
    gess_res <- gess_cmap(qSig=qsig_cmap, chunk_size=5000)
    score_col <- "scaled_score"
  }  
  if(gess_method == "LINCS"){
    query <- getDEGSig(drug, cell, Nup=Nup, Ndown=Ndown, refdb=refdb)
    qsig_lincs <- qSig(query = list(upset=query[[1]], downset=query[[2]]), 
                       gess_method="LINCS", refdb=refdb)
    gess_res <- gess_lincs(qsig_lincs, sortby="NCS", tau=tau)
    score_col <- "NCS"
  }
  if(gess_method == "gCMAP"){
    query <- getSig(drug, cell, refdb=refdb)
    qsig_gcmap <- qSig(query=query, gess_method="gCMAP", refdb=refdb)
    gess_res <- gess_gcmap(qsig_gcmap, higher=higher, lower=lower)
    score_col <- "effect"
  }
  if(gess_method == "Fisher"){
    query <- getSig(drug, cell, refdb=refdb)
    qsig_fisher <- qSig(query=query, gess_method="Fisher", refdb=refdb)
    gess_res <- gess_fisher(qSig=qsig_fisher, higher=higher, lower=lower)
    score_col <- "effect"
  }
  if(gess_method == "CORsub"){
    # draw expression value sig from lincs_expr
    query <- getSPsubSig(drug, cell, Nup=Nup, Ndown=Ndown)
    qsig_corsub <- qSig(query=query, gess_method="Cor", refdb="lincs_expr")
    gess_res <- gess_cor(qSig=qsig_corsub, method=method)
    score_col <- "cor_score"
  }
  if(gess_method == "CORall"){
    query <- getSig(drug, cell, refdb="lincs_expr")
    qsig_corall <- qSig(query=query, gess_method="Cor", refdb="lincs_expr")
    gess_res <- gess_cor(qSig=qsig_corall, method=method)
    score_col <- "cor_score"
  }
  gess_tb <- result(gess_res)
  drugs <- unique(gess_tb$pert)[seq_len(N_gess_drugs)]
  if(runFEA){
    if(fea_method == "dup_hyperG"){
      mf_res <- tsea_dup_hyperG(
        drugs=drugs, type="GO", ont="MF", 
        pvalueCutoff=pvalueCutoff, qvalueCutoff=qvalueCutoff,
        minGSSize=minGSSize, maxGSSize=maxGSSize)
      bp_res <- tsea_dup_hyperG(
        drugs=drugs, type="GO", ont="BP", 
        pvalueCutoff=pvalueCutoff, qvalueCutoff=qvalueCutoff,
        minGSSize=minGSSize, maxGSSize=maxGSSize)
      kegg_res <- tsea_dup_hyperG(
        drugs=drugs, type="KEGG",
        pvalueCutoff=pvalueCutoff, qvalueCutoff=qvalueCutoff,
        minGSSize=minGSSize, maxGSSize=maxGSSize)
    }
    if(fea_method == "mGSEA"){
      mf_res <- tsea_mGSEA(drugs=drugs, type="GO", ont="MF", exponent=1, 
                           nPerm=1000, pvalueCutoff=pvalueCutoff, minGSSize=minGSSize,
                           maxGSSize=maxGSSize)
      bp_res <- tsea_mGSEA(drugs=drugs, type="GO", ont="BP", 
                           pvalueCutoff=pvalueCutoff, minGSSize=minGSSize, maxGSSize=maxGSSize)
      kegg_res <- tsea_mGSEA(drugs=drugs, type="KEGG", pvalueCutoff=pvalueCutoff, 
                             minGSSize=minGSSize, maxGSSize=maxGSSize)
    } 
    if(fea_method == "mabs"){
      mf_res <- tsea_mabs(drugs=drugs, type="GO", ont="MF", 
                          pvalueCutoff=pvalueCutoff, minGSSize=minGSSize, maxGSSize=maxGSSize)
      bp_res <- tsea_mabs(drugs=drugs, type="GO", ont="BP", 
                          pvalueCutoff=pvalueCutoff, minGSSize=minGSSize, maxGSSize=maxGSSize)
      kegg_res <- tsea_mabs(drugs=drugs, type="KEGG", 
                            pvalueCutoff=pvalueCutoff, minGSSize=minGSSize, maxGSSize=maxGSSize)
    }
    if(fea_method == "hyperG"){
      mf_res <- dsea_hyperG(drugs=drugs, type="GO", ont="MF", 
                            pvalueCutoff=pvalueCutoff, qvalueCutoff=qvalueCutoff,
                            minGSSize=minGSSize, maxGSSize=maxGSSize)
      bp_res <- dsea_hyperG(drugs=drugs, type="GO", ont="BP", 
                            pvalueCutoff=pvalueCutoff, qvalueCutoff=qvalueCutoff,
                            minGSSize=minGSSize, maxGSSize=maxGSSize)
      kegg_res <- dsea_hyperG(drugs=drugs, type="KEGG", 
                              pvalueCutoff=pvalueCutoff, qvalueCutoff=qvalueCutoff,
                              minGSSize=minGSSize, maxGSSize=maxGSSize)
    }
    if(fea_method == "GSEA"){
      dl <- abs(gess_tb[[score_col]])
      names(dl) <- result(gess_res)$pert
      dl <- dl[dl>0]; dl <- dl[!duplicated(names(dl))]
      mf_res <- dsea_GSEA(drugList=dl, type="GO", ont="MF", exponent=1, 
                          nPerm=1000, pvalueCutoff=pvalueCutoff, minGSSize=minGSSize,
                          maxGSSize=maxGSSize)
      bp_res <- dsea_GSEA(drugList=dl, type="GO", ont="BP",
                          pvalueCutoff=pvalueCutoff, minGSSize=minGSSize, maxGSSize=maxGSSize)
      kegg_res <- dsea_GSEA(drugList=dl, type="KEGG", 
                            pvalueCutoff=pvalueCutoff, minGSSize=minGSSize, maxGSSize=maxGSSize)
    }
    mf_tb <- result(mf_res)
    bp_tb <- result(bp_res)
    kegg_tb <- result(kegg_res)
  }
  env_name <- paste0(env_dir, "/", drug, "_", cell, "_wf")
  res_dir <- paste0(env_name, "/results")
  if(!dir.exists(env_name)){
    dir.create(res_dir, recursive=TRUE)
  }
  write_tsv(gess_tb, paste0(res_dir, "/", gess_method, "_res.tsv"))
  
  if(runFEA){
    write_tsv(mf_tb, paste0(res_dir, "/", fea_method, "_mf_res.tsv"))
    write_tsv(bp_tb, paste0(res_dir, "/", fea_method, "_bp_res.tsv"))
    write_tsv(kegg_tb, paste0(res_dir, "/", fea_method, "_kegg_res.tsv"))
    
    if(GenerateReport){
      file.copy(system.file("extdata", "GESS_FEA_report.Rmd", package="signatureSearch"),
                paste0(env_name, "/GESS_FEA_report.Rmd"))
      Drug <- paste0(toupper(substr(drug, 1, 1)), substr(drug, 2, nchar(drug)))
      tx <- readLines(paste0(env_name, "/GESS_FEA_report.Rmd"))
      # replace tags in template rmd, such as drug, cell etc
      tx %<>% gsub("<Drug>", Drug, .) %>% gsub("<drug>", drug, .) %>%
        gsub("<cell>", cell, .) %>% gsub("<gess_method>", gess_method, .) %>%
        gsub("<fea_method>", fea_method, .) %>% gsub("<refdb>", refdb, .) %>%
        gsub("<N_gess_drugs>", N_gess_drugs, .)
      writeLines(tx, paste0(env_name, "/GESS_FEA_report.Rmd"))
      if(requireNamespace("rmarkdown", quietly=TRUE)){
        rmarkdown::render(paste0(env_name, "/GESS_FEA_report.Rmd"), quiet=TRUE)
      } }
  } else {
    if(GenerateReport){
      file.copy(system.file("extdata", "GESS_report.Rmd", package="signatureSearch"),
                paste0(env_name, "/GESS_report.Rmd"))
      Drug <- paste0(toupper(substr(drug, 1, 1)), substr(drug, 2, nchar(drug)))
      tx <- readLines(paste0(env_name, "/GESS_report.Rmd"))
      # replace tags in template rmd, such as drug, cell etc
      tx %<>% gsub("<Drug>", Drug, .) %>% gsub("<drug>", drug, .) %>%
        gsub("<cell>", cell, .) %>% gsub("<gess_method>", gess_method, .) %>%
        # gsub("<fea_method>", fea_method, .) %>% gsub("<refdb>", refdb, .) %>%
        gsub("<N_gess_drugs>", N_gess_drugs, .)
      writeLines(tx, paste0(env_name, "/GESS_report.Rmd"))
      if(requireNamespace("rmarkdown", quietly=TRUE)){
        rmarkdown::render(paste0(env_name, "/GESS_report.Rmd"), quiet=TRUE)
      } }
    
  }
  return(list(gess_tb=gess_tb, mf_tb=mf_tb, bp_tb=bp_tb, kegg_tb=kegg_tb))
}
