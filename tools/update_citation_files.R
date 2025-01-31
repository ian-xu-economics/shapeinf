update_citation_files <- function() {
  # Read the version number from DESCRIPTION
  desc <- read.dcf("DESCRIPTION")
  version <- desc[1, "Version"]
  
  # Update CITATION.cff
  cff_path <- "CITATION.cff"
  cff <- readLines(cff_path)
  cff <- sub("^version: .*", paste("version:", version), cff)
  writeLines(cff, cff_path)
  
  # Update inst/CITATION
  citation_path <- "inst/CITATION"
  citation <- readLines(citation_path)
  citation <- sub('version .*\\\"', paste0('version ', version, '"'), citation)
  writeLines(citation, citation_path)
  
  message("CITATION files updated to version ", version)
}

# Run the function
update_citation_files()
