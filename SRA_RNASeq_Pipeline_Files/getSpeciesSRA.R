#!/usr/bin/env Rscript

getSpeciesSRA <- function(SRR) {
  # # Bug testing
  # SRR <- "SRR1555623"
  # runInfo <- "../../proj_dir_2/Code/runInfo_Table.txt"
  
  runInfo <- "Code/runInfo_Table.txt"
  runInfoTable <- read.table(runInfo, sep = "\t", header = T, stringsAsFactors = F)
  if (! "ScientificName" %in% colnames(runInfoTable)) {
    return("unknown")
  }
  species <- runInfoTable$ScientificName[which(runInfoTable$Run == SRR)]
  if (species == "Homo sapiens") {
    species <- "human"
  } else if (species == "Mus musculus") {
    species <- "mouse"
  }
  
  return(species)
  
}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]

# Return result
res <- getSpeciesSRA(arg)
cat(res)


