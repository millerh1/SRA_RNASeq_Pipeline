#!/usr/bin/env Rscript

# Function to pull strandedness info from Salmon output file
getStrandednessInfo <- function(nameString) {
  js <- jsonlite::read_json(path = 
                              file.path("Results",
                                        "Salmon.out", 
                                        nameString, 
                                        "lib_format_counts.json"))
  libType <- js$expected_format
  res <- grep(libType, pattern = "S")
  if (length(res) != 0) {
    res2 <- grep(libType, pattern = "R")
    if (length(res2) != 0) {
      fin <- "reverse"
    } else {
      fin <- "yes"
    }
  } else {
    fin <- "no"
  }
  return(fin)
}

# Parse shell args
args <- commandArgs()
arg <- args[6]

# Return result
result <- getStrandednessInfo(arg)
cat(result)
