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
#' @param Signature The signature to perform signatureSearching on.
#' @param cellInfo A data table containing information about the cell types 
#' that make up the signature search database. This table must contain a column
#' labeled "cell_id" that contains matching cell type names in the signature 
#' search database.  Details about this file  can be found 
#' in the \code{cell_info2} table after load the `signatureSearch` package
#' and running `data("cell_info2")`
#' @param PertColName The name of the column from the signature search results 
#' table used to summarize results by.  
#' @param drug character(1) representing query drug name (e.g. vorinostat). 
#' This query drug should be included in the \code{refdb}
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
#' library(signatureSearch)
#' library(ExperimentHub); library(rhdf5)
#' eh <- ExperimentHub()
#' cmap <- eh[["EH3223"]]; cmap_expr <- eh[["EH3224"]]
#' lincs <- eh[["EH3226"]]; lincs_expr <- eh[["EH3227"]]
#' lincs2 <- eh[["EH7297"]]
#' h5ls(lincs2)
#' db_path <- system.file("extdata", "sample_db.h5", package = "signatureSearch")
#' library(SummarizedExperiment);
#' library(HDF5Array)
#' sample_db <- SummarizedExperiment(HDF5Array(db_path, name="assay"))
#' rownames(sample_db) <- HDF5Array(db_path, name="rownames")
#' colnames(sample_db) <- HDF5Array(db_path, name="colnames")
#' query_mat <- as.matrix(assay(sample_db[,"vorinostat__SKB__trt_cp"]))
#' query <- as.numeric(query_mat); names(query) <- rownames(query_mat)
#' upset <- head(names(query[order(-query)]), 150)
#' head(upset)
#' downset <- tail(names(query[order(-query)]), 150)
#' head(downset)
#' runWF(Signature = list(upset=upset, downset=downset),
#'       cellInfo = cellInfo2,
#'       PertColName = "pert_iname",
#'       drug = "vorinostat",
#'       refdb = lincs2,
#'       gess_method="LINCS",
#'       fea_method="dup_hyperG",
#'       N_gess_drugs=150,
#'       env_dir="./GESSWFResults",
#'       tau=FALSE,
#'       Nup=150,
#'       Ndown=150,
#'       higher=1,
#'       lower=-1,
#'       method="spearman",
#'       pvalueCutoff=1,
#'       qvalueCutoff=1,
#'       minGSSize=5,
#'       maxGSSize=500,
#'       runFEA=TRUE,
#'       GenerateReport= TRUE)
runWF <- function(Signature,cellInfo,PertColName = "pert_iname",drug,refdb,
                  gess_method="LINCS",fea_method="dup_hyperG",N_gess_drugs=150,env_dir=".",
                  tau=FALSE,Nup=150,Ndown=150,higher=1,lower=-1,
                  method="spearman",pvalueCutoff=1,qvalueCutoff=1,
                  minGSSize=5,maxGSSize=500,runFEA=TRUE,GenerateReport= TRUE){
  
  if(gess_method == "CMAP"){
    qsig_cmap <- qSig(query = Signature,
                      gess_method="CMAP", refdb=refdb)
    gess_res <- gess_cmap(qSig=qsig_cmap, chunk_size=5000)
    score_col <- "scaled_score"
  }
  if(gess_method == "LINCS"){
    qsig_lincs <- qSig(query = Signature,
                       gess_method="LINCS", refdb=refdb)
    gess_res <- gess_lincs(qsig_lincs, sortby="NCS", tau=tau)
    score_col <- "NCS"
  }
  if(gess_method == "gCMAP"){
    qsig_gcmap <- qSig(query=Signature, gess_method="gCMAP", refdb=refdb)
    gess_res <- gess_gcmap(qsig_gcmap, higher=higher, lower=lower)
    score_col <- "effect"
  }
  if(gess_method == "Fisher"){
    qsig_fisher <- qSig(query=Signature, gess_method="Fisher", refdb=refdb)
    gess_res <- gess_fisher(qSig=qsig_fisher, higher=higher, lower=lower)
    score_col <- "effect"
  }
  if(gess_method == "CORsub"){
    qsig_corsub <- qSig(query=Signature, gess_method="Cor", refdb="lincs_expr")
    gess_res <- gess_cor(qSig=qsig_corsub, method=method)
    score_col <- "cor_score"
  }
  if(gess_method == "CORall"){
    qsig_corall <- qSig(query=Signature, gess_method="Cor", refdb="lincs_expr")
    gess_res <- gess_cor(qSig=qsig_corall, method=method)
    score_col <- "cor_score"
  }
  gess_tb <- as.data.table(result(gess_res))
  
  #### Annotate GESS results ####
  sedb <- LINCSseLoad(DBpath = refdb)
  lincs_sig_info2 <- LINCSSigInfoGen(LINCSSummExp = sedb)
  data("lincs_pert_info2")
  AnnotSigInfo2 <- unique(merge(lincs_sig_info2, lincs_pert_info2[,1:2], by.x = "pert", by.y = "pert_id"))
  setnames(gess_tb, "pert", "pert_id", skip_absent=TRUE)
  gess_tb <- as.data.table(merge(gess_tb, lincs_pert_info2, by = "pert_id", all.x = TRUE))
  gess_tb <- gess_tb[order(gess_tb[[score_col]], decreasing = TRUE),]

  drugs <- unique(gess_tb[[PertColName]])[seq_len(N_gess_drugs)]
  if(runFEA){
    print("Performing FEA analysis")
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
  
  #### Calculate drug ranking by cell type ####
  print("Calculating drug ranks by cell type")
  #### Classify ranks across cell types ####
  if(gess_method == "LINCS"){
    CellCat <- GESSAttributeCatalog(ClasifyDT = gess_tb, RowFeature = "pert_iname", ColFeature = "cell",
                                    ValueCol = "Rank", method = "mean", addOrderingRow = FALSE,
                                    Rescore = TRUE, ScoreCol = score_col, cellInfo = cellInfo)
  }
  
  #### Calculate expression level for each signature gene across all results ####
  print("Obtaining Signature Expression data")
  degMat <- getAllSig(refdb=refdb, gess_tb=gess_tb, Signature=Signature, method="mean")
  
  print("Writing report")
  env_name <- paste0(env_dir)
  res_dir <- paste0(env_name, "/results")
  if(!dir.exists(env_name)){ dir.create(res_dir, recursive=TRUE) }
  fwrite(gess_tb, paste0(res_dir, "/", gess_method, "_res.xls"), row.names=FALSE, quote=FALSE, sep="\t")
  fwrite(CellCat, paste0(res_dir, "/", "ResultRankByCell.xls"), row.names=FALSE, quote=FALSE, sep="\t")
  fwrite(degMat, paste0(res_dir, "/", "GESExpressionLevel.xls"), row.names=FALSE, quote=FALSE, sep="\t")
  
  if(runFEA){
    fwrite(mf_tb, paste0(res_dir, "/", fea_method, "_mf_res.xls"), row.names=FALSE, quote=FALSE, sep="\t")
    fwrite(bp_tb, paste0(res_dir, "/", fea_method, "_bp_res.xls"), row.names=FALSE, quote=FALSE, sep="\t")
    fwrite(kegg_tb, paste0(res_dir, "/", fea_method, "_kegg_res.xls"), row.names=FALSE, quote=FALSE, sep="\t")
    
    if(GenerateReport){
      file.copy(system.file("extdata", "GESS_FEA_report.Rmd", package="signatureSearch"),
                paste0(env_name, "/GESS_FEA_report.Rmd"))
      Drug <- paste0(toupper(substr(drug, 1, 1)), substr(drug, 2, nchar(drug)))
      tx <- readLines(paste0(env_name, "/GESS_FEA_report.Rmd"))
      # replace tags in template rmd, such as drug, cell etc
      tx %<>% gsub("<Drug>", Drug, .) %>% gsub("<drug>", drug, .) %>%
        gsub("<gess_method>", gess_method, .) %>%
        gsub("<fea_method>", fea_method, .) %>% gsub("<refdb>", refdb, .) %>%
        gsub("<N_gess_drugs>", N_gess_drugs, .) %>%
        gsub("<env_name>", env_name, .)
      writeLines(tx, paste0(env_name, "/GESS_FEA_report.Rmd"))
      if(requireNamespace("rmarkdown", quietly=TRUE)){
        rmarkdown::render(paste0(env_name, "/GESS_FEA_report.Rmd"), quiet=TRUE)
      }
    }
  } else {
    if(GenerateReport){
      file.copy(system.file("extdata", "GESS_report.Rmd", package="signatureSearch"),
                paste0(env_name, "/GESS_report.Rmd"))
      Drug <- paste0(toupper(substr(drug, 1, 1)), substr(drug, 2, nchar(drug)))
      tx <- readLines(paste0(env_name, "/GESS_report.Rmd"))
      # replace tags in template rmd, such as drug, cell etc
      tx %<>% gsub("<Drug>", Drug, .) %>% gsub("<drug>", drug, .) %>%
        gsub("<gess_method>", gess_method, .) %>%
        gsub("<refdb>", refdb, .) %>%
        gsub("<N_gess_drugs>", N_gess_drugs, .) %>%
        gsub("<env_name>", env_name, .)
      writeLines(tx, paste0(env_name, "/GESS_report.Rmd"))
      if(requireNamespace("rmarkdown", quietly=TRUE)){
        rmarkdown::render(paste0(env_name, "/GESS_report.Rmd"), quiet=TRUE)
      } }
  }
  print("Analysis completed")
  if(runFEA){
    return(list(gess_tb=gess_tb, CellGESS=CellCat, DEG=degMat, mf_tb=mf_tb, bp_tb=bp_tb, kegg_tb=kegg_tb))
  } else {
    return(list(gess_tb=gess_tb, CellGESS=CellCat, DEG=degMat))
  }
}

#' @import SummarizedExperiment
#' @import HDF5Array
LINCSseLoad <- function(DBpath){
  sedb <- SummarizedExperiment(HDF5Array(DBpath, name="assay"))
  rownames(sedb) <- HDF5Array(DBpath, name="rownames")
  colnames(sedb) <- HDF5Array(DBpath, name="colnames")
  return(sedb)}

#' @importFrom stringr str_split
#' @import SummarizedExperiment
LINCSSigInfoGen <- function(LINCSSummExp = sedb){
  spl <- str_split(colnames(LINCSSummExp), "__")
  pert <- NULL; for(i in 1:length(spl)){pert[i] <- spl[[i]][1]}
  cell <- NULL; for(i in 1:length(spl)){cell[i] <- spl[[i]][2]}
  pert_type <- NULL; for(i in 1:length(spl)){pert_type[i] <- spl[[i]][3]}
  lincs_sig_info2 <- data.frame(pert = pert, cell = cell, pert_type =pert_type)
  return(lincs_sig_info2)}

#' @import data.table
GESSAttributeCatalog <- function(ClasifyDT, RowFeature, ColFeature, ValueCol, method, addOrderingRow,
                                 Rescore, ScoreCol = "NCS", cellInfo){
  #### summarize results by cell type ####
  gess_tb2 <- ClasifyDT[!ClasifyDT[[ScoreCol]] == 0,]
  gess_tb2 <- gess_tb2[order(gess_tb2[[ScoreCol]], decreasing = TRUE),][,Rank := 1:nrow(gess_tb2)][]
  #### Cross reference to cell types in GESS results ####
  cell_info2Used <- as.data.table(cellInfo[cellInfo$cell_id %in% unique(gess_tb2$cell), ]) %>% setnames("cell_id", "cell")
  #### Annotate with the representation of cell types in the GESS results ####
  cdt <- data.table(table(gess_tb2$cell)) %>% setnames(c("V1", "N"), c("cell", "count"))
  cell_info2Used <- merge(cell_info2Used, cdt, by = "cell")
  gess_tb2 <- merge(gess_tb2, cell_info2Used, by = "cell") # dim(gess_tb2); dim(ClassificationDT)
  gess_tb2 <- gess_tb2[order(Rank, decreasing = FALSE),]
  #### Re-score by attribute ####
  if(Rescore){
    tempDT <- data.table()
    feat <- unique(gess_tb2[[ColFeature]])
    feat <- feat[!feat == ""]
    for(c in 1:length(feat)){
      t <- gess_tb2[gess_tb2[[ColFeature]] == feat[c],]
      t <- t[,.SD[which.min(Rank)], by = RowFeature]
      t <- t[order(t[[ScoreCol]], decreasing = TRUE),]
      t[,Rank := 1:nrow(t)]
      tempDT <- rbind(tempDT, t)
    }
    gess_tb2 <- tempDT
  }
  #### Set up matrix ####
  Rows <- unique(gess_tb2[[RowFeature]])
  Rows <- Rows[!Rows == ""]
  Rows <- Rows[!is.na(Rows)]
  Cols <- unique(gess_tb2[[ColFeature]])
  Cols <- Cols[!Cols == ""]
  Cols <- Cols[!is.na(Cols)]
  mat <- matrix(NA, nrow = length(Rows), ncol = length(Cols))
  rownames(mat) <- Rows
  colnames(mat) <- Cols
  #### Compile values ####
  for(a in 1:length(Rows)){
    temp <- gess_tb2[gess_tb2[[RowFeature]] == Rows[a], ]
    for(b in 1:length(Cols)){
      if(method == "mean"){ mat[Rows[a], Cols[b]] <- mean(temp[temp[[ColFeature]] == Cols[b],][[ValueCol]]) }
      if(method == "median"){ mat[Rows[a], Cols[b]] <- median(temp[temp[[ColFeature]] == Cols[b],][[ValueCol]]) }
      if(method == "sum"){ mat[Rows[a], Cols[b]] <- sum(temp[temp[[ColFeature]] == Cols[b],][[ValueCol]]) }
    } }
  logic <- !is.na(mat)
  #### order matrix ####
  rowSums <- apply(logic, 1, sum)
  ColSums <- apply(logic, 2, sum)
  df <- as.data.frame(mat)
  df$RowSums <- rowSums
  if(addOrderingRow){
    ColumnSums <- 0; names(ColumnSums) <- "ColSums"
    ColSums <- c(ColSums, ColumnSums)
    coladd <- t(data.frame(RowSums))
    df <- rbind(df, coladd)
  }
  df <- df[order(df$RowSums, decreasing = TRUE),  ]
  df$Pert_iname <- rownames(df)
  df <- as.data.table(df[,c(ncol(df), c(1:(ncol(df) -1)))])
  return(df)
}

#' @import org.Hs.eg.db
#' @import AnnotationHub
#' @import SummarizedExperiment
#' @import HDF5Array
#' @import data.table
getAllSig <- function(refdb, gess_tb, Signature, method){
  fullSig <- as.character(unlist(Signature))
  if (is.character(refdb)) {
    refdb <- signatureSearch:::determine_refdb(refdb)
    refse <- SummarizedExperiment(HDF5Array(refdb, name = "assay"))
    rownames(refse) <- HDF5Array(refdb, name = "rownames")
    colnames(refse) <- gsub("__trt_cp", "", HDF5Array(refdb, name = "colnames"))
    #### Obtain expression values ####
    index <- unique(paste(gess_tb$pert_iname))
    expDT <- rbindlist(lapply(seq_along(index), function(x){
      temp <- gess_tb[pert_iname == index[x],]
      ids <- unique(paste(temp$pert_id, temp$cell, sep ="__"))
      trt2 <- intersect(ids, colnames(refse))
      cmp_mat <- as.matrix(assay(refse[, trt2]))
      if(ncol(cmp_mat) > 1){
        matSub <- cmp_mat[row.names(cmp_mat) %in% fullSig,]
        if(method == "mean"){ dt <- as.data.table(t(data.frame(rowMeans(matSub))))  }
      } else {dt <- as.data.table(t(data.frame(cmp_mat[row.names(cmp_mat) %in% fullSig,])) )   }
      return(dt)  }) )
    #### Map ENTREZID to Gene Names ####
    geneOnt <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys=colnames(expDT), columns=c("SYMBOL"), keytype="ENTREZID")) # , "ENSEMBL"
    geneOnt <- geneOnt[match(ENTREZID, colnames(expDT)),]
    geneOnt[is.na(SYMBOL),]$SYMBOL <- "NoName"
    setnames(expDT, as.character(geneOnt$ENTREZID), geneOnt$SYMBOL)
    #### Add drug names ####
    cmp_DT <- expDT[, pert := index][,c(ncol(expDT), 1:ncol(expDT)-1), with = FALSE]
    return(cmp_DT)
  } else {
    message(paste("Please set refdb as one of", "'lincs', 'lincs_expr', 'cmap' or 'cmap_expr', "),
            "or path to an HDF5 file representing reference database!")
  }
}


