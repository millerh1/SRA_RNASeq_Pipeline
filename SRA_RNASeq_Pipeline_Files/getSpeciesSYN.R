#!/usr/bin/env Rscript

getSpeciesSYN <- function(line, runInfo) {
  
  # # # Bug testing
  # runInfo <- "~/TESTGLIA/Code/runInfo_Table_2019-07-22_11.52.06.txt"
  
  mapFile <- runInfo
  map <- read.table(mapFile, stringsAsFactors = F, sep = "\t", header = T)
  species <- map$species[which(map$Run == line)]
  species <- tolower(species)
  format <- unique(map$fileFormat)
  return(species)
}


# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]
arg2 <- args[2]
result <- getSpeciesSYN(arg, arg2)
cat(result)
