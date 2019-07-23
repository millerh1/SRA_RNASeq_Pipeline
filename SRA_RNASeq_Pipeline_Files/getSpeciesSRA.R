#!/usr/bin/env Rscript

getSpeciesSRA <- function(SRR, runInfo) {
  # # # Bug testing
  # SRR <- "SRR5527012"
  # runInfo <- "../../Bishop.lab/EWS_CTR_All_Cells/Code/runInfo_Table_2019-07-23_09.07.29.txt"
  # system(paste0("dos2unix ", runInfo))
  # cat("\n", file = runInfo, append = TRUE)
  
  runInfoTable <- read.csv(runInfo, header = T, stringsAsFactors = F)
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
arg2 <- args[2]

# Return result
res <- getSpeciesSRA(arg, arg2)
cat(res)


