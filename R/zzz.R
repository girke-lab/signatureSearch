#' @import ExperimentHub
eh <- NULL
GO_DATA <- NULL
GO_DATA_drug <- NULL
.onLoad <- function(libname, pkgname){
    # .some_internal_global_variable <<- 55L
    eh <<- tryCatch({
        cache <- tools::R_user_dir("ExperimentHub", which="cache")
        packageStartupMessage("The ExperimentHub cache is at ", cache)
        # setExperimentHubOption("CACHE", cache)
        suppressMessages(ExperimentHub())
        }, error=function(e){
            stop(
        "The ExperimentHub cache is corrupt, please refer to 
  https://bioconductor.org/packages/release/bioc/vignettes/AnnotationHub/inst/doc/TroubleshootingTheCache.html for debugging.")
            })
    # GO_DATA <- get_GO_data(OrgDb, ont, keytype="SYMBOL")
    # download GO_DATA.rds from AnnotationHub to save time by avoiding 
    # building GO_DATA from scratch
    GO_DATA <<- validLoad("EH3231")
    
    # GO_DATA_drug <- get_GO_data_drug(OrgDb = "org.Hs.eg.db", 
    #                                  ont, keytype="SYMBOL")
    # download GO_DATA_drug.rds 
    GO_DATA_drug <<- validLoad("EH3232")
}

#' @importFrom BiocGenerics fileName
validLoad <- function(ehid){
    tryCatch(suppressMessages(eh[[ehid]]), 
             error=function(e){
                 unlink(fileName(eh[ehid]))
                 eh[[ehid]]})
}