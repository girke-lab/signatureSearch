## Function to get PubChem CIDs by name
name2pubchemCID <- function(x) {
  #print(x)
  url <- paste0('"https://pubchem.ncbi.nlm.nih.gov/compound/', x, '"')
  system2("wget", paste(url, '-O tmp_file --quiet'))
  pc_page <- readLines("tmp_file", warn=FALSE)
  unlink("tmp_file")
  cid_line <- pc_page[grepl("meta name=\"pubchem_uid_value\" content=", pc_page)][1]
  cid <- gsub("(^.*content=\")|\">$", "", cid_line)
  if(is.na(cid[1])){
    return("Not Found")
  } else {
    return(cid)
  }
}
## Usage:
# ids <- c("rhodomyrtoxin-b", "SD-169", "not_avail_test", "baccatin-III")
# vapply(ids, function(i) name2pubchemCID(i), character(1))

getPertInfo <- function(pert_iname, user_key){
  pert_info = NULL
  col = c("canonical_smiles", "description", "inchi_key", "inchi_string", "moa", 
          "molecular_formula", "pert_id", "pert_iname", "pert_type", "pubchem_cid", "target")

  for(iname in pert_iname){
    json <- system(paste0('curl -X GET --globoff --header "Accept: application/json" --header "user_key: ',
    user_key, '" "https://api.clue.io/api/perts?filter=%7B%22where%22%3A%7B%22pert_iname%22%3A%22',
    iname,'%22%7D%7D"'), intern = TRUE, ignore.stderr = TRUE)
    
    if(length(json)==0){
      message("Curl error happens on ", iname)
      next()
    }
    if(grepl("error",json) | json=="[]"){
      message("API error happens on ", iname)
      next()
    }
    require(jsonlite)
    jsonr <- fromJSON(json)
    # print(i)
    # i=i+1
    # rownames(jsonr) <- pert_id
    for(item in col){
      if(! item %in% colnames(jsonr)){
        new <- data.frame(NA)
        colnames(new) <- item
        jsonr = cbind(jsonr, new)
      }
    }
    pert_info <- rbind(pert_info, jsonr[,col])
  }
  return(pert_info)
}
