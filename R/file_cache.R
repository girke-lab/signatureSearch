#' @import BiocFileCache
#' @importFrom rappdirs user_cache_dir
.get_cache <- function(){
    cache <- rappdirs::user_cache_dir(appname="signatureSearch")
    BiocFileCache::BiocFileCache(cache)
}

download_data_file <- function(url, rname, verbose=FALSE){
    fileURL <- url
    bfc <- suppressMessages(.get_cache())
    rid <- bfcquery(bfc, rname, "rname", exact=TRUE)$rid
    if (!length(rid)) {
      if(verbose)
        message(paste("Downloading",rname, "from", url))
      rid <- names(bfcadd(bfc, rname, fileURL))
    }
    if (!isFALSE(bfcneedsupdate(bfc, rid)))
      bfcdownload(bfc, rid)
    bfcrpath(bfc, rids = rid)
}
