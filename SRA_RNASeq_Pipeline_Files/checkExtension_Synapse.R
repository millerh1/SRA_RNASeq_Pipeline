#!/usr/bin/env Rscript

checkExtension_Synapse <- function(runInfo) {
  
  # # # Bug testing
  # runInfo <- "~/TESTGLIA/Code/runInfo_Table_2019-07-22_11.52.06.txt"
  
  mapFile <- runInfo
  map <- read.table(mapFile, stringsAsFactors = F, sep = "\t", header = T)
  format <- unique(map$fileFormat)
  if (length(format == 1)) {
    res <- format
  } else {
    res <- "ERROR: cannot currently handle mixed file types. Please select accessions with only one format type (fastq or bam)."
  }
  return(res)
}


# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]
result <- checkExtension_Synapse(arg)
cat(result)
