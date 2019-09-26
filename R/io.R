########################################
### GCT class and method definitions ###
########################################

#' An S4 Class to Represent a GCT Object
#' @name GCT object
#' @slot mat a numeric matrix
#' @slot rid a character vector of row ids
#' @slot cid a character vector of column ids
#' @slot rdesc a \code{data.frame} of row descriptors
#' @slot rdesc a \code{data.frame} of column descriptors
#' @slot src a character indicating the source (usually file path) 
#' of the data
#' @description The GCT class serves to represent annotated
#'   matrices. The \code{mat} slot contains the numeric matrix data and the
#'   \code{rdesc} and \code{cdesc} slots contain data frames with
#'   annotations about the rows and columns, respectively
#'   
#' @seealso \code{\link{parse_gctx}}
methods::setClass("GCT",
         methods::representation(
             mat = "matrix",
             rid = "character",
             cid = "character",
             rdesc = "data.frame",
             cdesc = "data.frame",
             version = "character",
             src = "character"
         )
)


# set up methods for checking GCT validity
methods::setValidity("GCT",
  function(object) {
    # check whether dimensions of various
    # slots are in sync
    nrows <- nrow(object@mat)
    ncols <- ncol(object@mat)
    if (nrows != length(object@rid)) {
      return("rid must be the same length as number of matrix rows")
    }
    if (ncols != length(object@cid)) {
      return("cid must be the same length as number of matrix columns")
    }
    if (length(object@cid) > length(unique(object@cid))) {
      return("cid must be unique")
    }
    if (length(object@rid) > length(unique(object@rid))) {
      return("rid must be unique")
    }
    if (nrow(object@cdesc) != ncols & nrow(object@cdesc) != 0) {
      return("cdesc must either have 0 rows or the same number of rows as 
             matrix has columns")
    }
    if (nrow(object@rdesc) != nrows & nrow(object@rdesc) != 0) {
      return("rdesc must either have 0 rows or the same number of rows as 
             matrix has rows")
    }
    else {
      return(TRUE)
    }
  }
)


#### define some helper methods for parsing gctx files ###

# Adjust the data types for columns of a meta data frame
fix.datatypes <- function(meta) {
    for (field.name in names(meta)) {
        # get the field values
        field <- meta[[field.name]]
        # check if it's numeric. data may come in as a string
        # but actually contains numeric values. if so, as.numeric
        # will not result in a vector of NA values
        field.as.numeric <- suppressWarnings(as.numeric(field))
        if (!any(is.na(field.as.numeric))) {
          field <- field.as.numeric
        }
        if (is.numeric(field)) {
            # check if it's an integer. data may be floats but
            # if we coerce to an integer and the difference from
            # original values is zero, that means data are actually
            # integers. integer conversion will return NA if there
            # are any issues.
            field.as.integer <- suppressWarnings(as.integer(field))
            if (!any(is.na(field.as.integer))) {
              # integer conversion was fine, lets see if the
              # values are altered
              diffs <- field - field.as.integer
              if (all(diffs == 0)) {
                # converting to integer didn't change value,
                # set field to integer values
                field <- field.as.integer
              }
            }
        }
        # insert back into the annotations
        meta[[field.name]] <- field
    }
    return(meta)
}


# Parse row or column metadata from GCTX files
read.gctx.meta <- function(gctx_path, dimension="row", ids=NULL, 
                           set_annot_rownames=TRUE) {
  if (!file.exists(gctx_path)) {
    stop(paste(gctx_path, "does not exist"))
  }
  if (dimension=="column") dimension <- "col"
  if (!(dimension %in% c("row", "col"))) {
    stop("dimension can be either row or col")
  }
  if (dimension == "row") {
    name <- "0/META/ROW"
  } else {
    name <- "0/META/COL"
  }
  raw_annots <- rhdf5::h5read(gctx_path, name=name) # returns a list
  fields <- names(raw_annots)
  # define an empty data frame of the correct dimensions
  annots <-  data.frame(matrix(nrow=length(raw_annots[[fields[1]]]), 
                               ncol=length(fields)))
  names(annots) <- fields
  # loop through each field and fill the annots data.frame
  for (i in seq_along(fields)) {
    field <- fields[i]
    # remove any trailing spaces
    # and cast as vector
    annots[,i] <- as.vector(gsub("\\s*$", "", raw_annots[[field]], perl=TRUE))
  } 
  annots <- fix.datatypes(annots)
  # subset to the provided set of ids, if given
  if (is.null(ids)) {
    ids <- as.character(annots$id)
  } else {
    ids <- ids
  }
  # make sure annots row ordering matches that of ids
  annots <- subset_to_ids(annots, ids)
  annots$id <- as.character(annots$id)
  # use the id field to set the rownames
  if (set_annot_rownames) {
    rownames(annots) <- annots$id
  }
  return(annots)
}


# Read GCTX row or column ids
read.gctx.ids <- function(gctx_path, dimension="row") {
  if (!file.exists(gctx_path)) {
    stop(paste(gctx_path, "does not exist"))
  }
  if (dimension=="column") dimension <- "col"
  if (!(dimension %in% c("row", "col"))) {
    stop("dimension can be either row or col")
  }
  if (dimension == "row") {
    name <- "0/META/ROW/id"
  } else {
    name <- "0/META/COL/id"
  }
  # remove any spaces
  ids <- gsub("\\s*$", "", rhdf5::h5read(gctx_path, name=name), perl=TRUE)
  # cast as character
  ids <- as.character(ids)
  return(ids)
}

# Return a subset of requested GCTX row/colum ids
# out of the universe of all ids
process_ids <- function(ids, all_ids, type="rid") {
  if (!is.null(ids)) {
    if (is.numeric(ids)) {
      # is it numeric?
      idx <- ids
      is_invalid_idx <- (idx > length(all_ids)) | (idx <= 0)
      invalid_idx <- idx[is_invalid_idx]
      if (all(is_invalid_idx)) {
        stop(paste("none of the requested", type, 
                   "indices were found in the dataset"))
      }
      if (any(is_invalid_idx)) {
        # requested indices are outside of the possible range
        warning(paste("the following ", type, 
                      " were are outside possible range and will be ignored:\n",
                      paste(invalid_idx, collapse="\n"), sep="")) 
      }
      idx <- idx[!is_invalid_idx]
    } else {
      # assume its a character
      idx <- match(ids, all_ids)
      if (all(is.na(idx))) {
        stop(paste("none of the requested", type, "were found in the dataset"))
      }
      if (any(is.na(idx))) {
        ids_not_found <- ids[is.na(idx)]
        warning(paste("the following ", type, 
                      " were not found and will be ignored:\n",
                      paste(ids_not_found, collapse="\n"), sep=""))
      }
      idx <- idx[!is.na(idx)]
    }
  } else {
    # ids were null, just return an index vector
    # allong all_ids
    idx <- seq_along(all_ids)
  }
  # subset the character ids to the ones we want
  id_keep <- as.character(all_ids[idx])
  return(list(idx=idx, ids=id_keep))
}


# define the initialization method for the GCT class
methods::setMethod("initialize",
          signature = "GCT",
          definition = function(.Object, mat=NULL, rdesc=NULL, cdesc=NULL, 
                                src=NULL, rid=NULL, cid=NULL, 
                                set_annot_rownames=FALSE,
                                matrix_only=FALSE) {
            # if we were supplied a matrix and annotations, use them
            if (!is.null(mat)) {
              .Object@mat <- mat
              # if given rid and cid, use those as well
              if (!is.null(rid)) {
                .Object@rid <- rid
              } else {
                .Object@rid <- rownames(mat)
              }
              if (!is.null(cid)) {
                .Object@cid <- cid
              } else {
                .Object@cid <- colnames(mat)
              }
            }
            if (!is.null(rdesc)) {
              .Object@rdesc <- rdesc
            }
            if (!is.null(cdesc)) {
              .Object@cdesc <- cdesc
            } else if (!is.null(src)) {
              # we were not given a matrix, were we given a src file?
              # check to make sure it's either .gct or .gctx
              if (! (grepl(".gct$", src) || grepl(".gctx$", src) ))
                stop("Either a .gct or .gctx file must be given")
              if (grepl(".gct$", src)) {
                if ( ! is.null(rid) || !is.null(cid) )
                  warning(paste("rid and cid values may only be given for 
                                .gctx files, not .gct files\n",
                                "ignoring"))
                # parse the .gct
                .Object@src = src
                # get the .gct version by reading first line
                .Object@version = scan(src, what = "", nlines = 1, sep = "\t", 
                                       quiet = TRUE)[1]
                # get matrix dimensions by reading second line
                dimensions = scan(src, what = double(0), nlines = 1, skip = 1, 
                                  sep = "\t", quiet = TRUE)
                nrmat = dimensions[1]
                ncmat = dimensions[2]
                if (length(dimensions)==4) {
                  # a #1.3 file
                  message("parsing as GCT v1.3")
                  nrhd <- dimensions[3]
                  nchd <- dimensions[4]
                } else {
                  # a #1.2 file
                  message("parsing as GCT v1.2")
                  nrhd <- 0
                  nchd <- 0
                }
                message(paste(src, nrmat, "rows,", ncmat, "cols,", nrhd, 
                              "row descriptors,", nchd, "col descriptors"))
                # read in header line
                header = scan(src, what = "", nlines = 1, skip = 2, sep = "\t", 
                              quote = NULL, quiet = TRUE)
                # construct row header and column id's from the header line
                if ( nrhd > 0 ) {
                  rhd <- header[2:(nrhd+1)]
                  cid <- header[-(nrhd+1):-1]
                  col_offset <- 1
                }
                else {
                  if (any(grepl("description", header, ignore.case=T))) {
                    # check for presence of description column in v1.2 files
                    col_offset <- 2
                  } else {
                    col_offset <- col_offset <- 1
                  }
                  rhd = NULL
                  cid = header[(1+col_offset):length(header)]
                }
                # read in the next set of headers (column annotations) and 
                # shape into a matrix
                if ( nchd > 0 ) {
                  header = scan(src, what = "", nlines = nchd, skip = 3, 
                                sep = "\t", quote = NULL, quiet = TRUE)
                  header = matrix(header, nrow = nchd, 
                                  ncol = ncmat + nrhd + 1, byrow = TRUE)
                  # extract the column header and column descriptions
                  chd = header[,1]
                  cdesc = header[,-(nrhd+1):-1]
                  # need to transpose in the case where there's only one column 
                  # annotation
                  if ( nchd == 1 )
                    cdesc = t(cdesc)
                }
                else {
                  chd = NULL
                  cdesc = data.frame()
                }
                # read in the data matrix and row descriptions, shape into a
                #  matrix
                mat = scan(src, what = "", nlines = nrmat, skip = 3 + nchd, 
                           sep = "\t", quote = NULL, quiet = TRUE)
                mat = matrix(mat, nrow = nrmat, 
                             ncol = ncmat + nrhd + col_offset, 
                             byrow = TRUE)
                # message(paste(dim(mat), collapse="\t"))
                # Extract the row id's row descriptions, and the data matrix
                rid = mat[,1]
                if ( nrhd > 0 ) {
                  # need as.matrix for the case where there's only one row 
                  # annotation
                  rdesc = as.matrix(mat[,2:(nrhd + 1)])
                  mat = matrix(as.numeric(mat[,-(nrhd + 1):-1]),
                               nrow = nrmat, ncol = ncmat)
                }
                else {
                  rdesc = data.frame()
                  mat = matrix(as.numeric(mat[, (1+col_offset):ncol(mat)]), 
                               nrow = nrmat, ncol = ncmat)
                }
                # assign names to the data matrix and the row and column 
                # descriptions
                # message(paste(dim(mat), collapse="\t"))
                dimnames(mat) = list(rid, cid)
                if ( nrhd > 0 ) {
                  dimnames(rdesc) = list(rid,rhd)
                  rdesc = as.data.frame(rdesc, stringsAsFactors = FALSE)
                }
                if ( nchd > 0 ) {
                  cdesc = t(cdesc)
                  dimnames(cdesc) = list(cid,chd)
                  cdesc = as.data.frame(cdesc, stringsAsFactors = FALSE)
                }
                # assign to the GCT slots
                .Object@mat = mat
                .Object@rid = rownames(mat)
                .Object@cid = colnames(mat)
                if (!matrix_only) {
                  # return annotations as well as matrix
                  .Object@rdesc = fix.datatypes(rdesc)
                  .Object@cdesc = fix.datatypes(cdesc)
                  # add id columns to rdesc and cdesc
                  .Object@rdesc$id <- rownames(.Object@rdesc)
                  .Object@cdesc$id <- rownames(.Object@cdesc)
                }
              }
              else { 
                # parse the .gctx
                message(paste("reading", src))
                .Object@src = src
                # get all the row and column ids
                all_rid <- read.gctx.ids(src, dimension="row")
                all_cid <- read.gctx.ids(src, dimension="col")
                # if rid or cid specified, read only those rows/columns
                # if already numeric, use as is
                # else convert to numeric indices
                processed_rids <- process_ids(rid, all_rid, type="rid")
                processed_cids <- process_ids(cid, all_cid, type="cid")
                # read the data matrix
                .Object@mat <- rhdf5::h5read(src, name="0/DATA/0/matrix",
                                      index=list(processed_rids$idx,
                                                 processed_cids$idx))
                # set the row and column ids, casting as characters
                .Object@rid <- processed_rids$ids
                .Object@cid <- processed_cids$ids
                rownames(.Object@mat) <- processed_rids$ids
                colnames(.Object@mat) <- processed_cids$ids
                # get the meta data
                if (!matrix_only) {
                  .Object@rdesc <- read.gctx.meta(src, dimension="row", 
                                                  ids=processed_rids$ids,
                                    set_annot_rownames=set_annot_rownames)
                  .Object@cdesc <- read.gctx.meta(src, dimension="col", 
                                                  ids=processed_cids$ids,
                                     set_annot_rownames=set_annot_rownames)
                }
                else {
                  .Object@rdesc <- data.frame(id=.Object@rid, 
                                              stringsAsFactors = FALSE)
                  .Object@cdesc <- data.frame(id=.Object@cid, 
                                              stringsAsFactors = FALSE)
                }
                # close any open handles and return the object
                if(utils::packageVersion('rhdf5') < "2.23.0") {
                    rhdf5::H5close()
                } else {
                  rhdf5::h5closeAll()
                }
                message("done")
              }
            }
            # finally, make sure object is valid before returning
            ok <- methods::validObject(.Object)
            return(.Object)
          }
)


#' Parse a GCTX file into the R workspace as a GCT object
#' @title Parse GCTX
#' @param fname character(1), path to the GCTX file on disk
#' @param rid either a vector of character or integer
#'   row indices or a path to a grp file containing character
#'   row indices. Only these indices will be parsed from the
#'   file.
#' @param cid either a vector of character or integer
#'   column indices or a path to a grp file containing character
#'   column indices. Only these indices will be parsed from the
#'   file.
#' @param set_annot_rownames boolean indicating whether to set the
#'   rownames on the row/column metadata data.frames. Set this to 
#'   false if the GCTX file has duplicate row/column ids.
#' @param matrix_only boolean indicating whether to parse only
#'   the matrix (ignoring row and column annotations)
#' @return gct object
#' @family GCTX parsing functions
#' @examples 
#' gctx <- system.file("extdata", "test_sample_n2x12328.gctx", 
#'                     package="signatureSearch")
#' gct <- parse_gctx(gctx)
#' @export
parse_gctx <- function(fname, rid=NULL, cid=NULL, set_annot_rownames=FALSE, 
                       matrix_only=FALSE) {
    ds <- methods::new("GCT",
              src = fname,
              rid = rid,
              cid = cid,
              set_annot_rownames = set_annot_rownames,
              matrix_only = matrix_only)
    return(ds)
}

# Do a robust \code{\link{data.frame}} subset to a set of ids
subset_to_ids <- function(df, ids) {
  # helper function to do a robust df subset
  check_colnames("id", df)
  newdf <- data.frame(df[match(ids, df$id), ])
  names(newdf) <- names(df)
  return(newdf)
}

# Check whether \code{test_names} are columns in the \code{\link{data.frame}} df
check_colnames <- function(test_names, df, throw_error=TRUE) {
  # check whether test_names are valid names in df
  # throw error if specified
  diffs <- setdiff(test_names, names(df))
  if (length(diffs) > 0) {
    if (throw_error) {
      stop(paste("the following column names are not found in", 
                 deparse(substitute(df)), ":",
                 paste(diffs, collapse=" "), "\n"))
    } else {
      return(FALSE)
    }
  } else {
    return(TRUE)
  }
}
