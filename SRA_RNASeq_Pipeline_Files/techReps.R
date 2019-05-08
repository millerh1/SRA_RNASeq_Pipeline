#!/usr/bin/env Rscript

techReps <- function(SRA, runInfoFile) {
  runInfo <- read.table(runInfoFile, sep = "\t", fill = T)
  if (SRA %in% runInfo$Run) {
    SRX <- runInfo$Experiment[which(runInfo$Run == SRA)]
    SRXdups <- runInfo$Experiment[which(duplicated(runInfo$Experiment))]
    if (SRX %in% SRXdups) {
      res <- 'yes'
      SRAs <- runInfo$Run[which(runInfo$Experiment == SRX)]
      ind <- which(SRAs == SRA)
      if (ind == 1) {
        SRAs <- as.matrix(SRAs)
        write.table(SRAs, file = "Data/tmp/sraTempTable.txt", quote = F, row.names = F)
      } else {
        res <- "skip"
      }
    } else {
      res <- 'no'
    }
  } else {
    res <- 'no'
  }
  return(res)
}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]
argTable <- args[2]

# Return result
result <- techReps(arg, argTable)
cat(result)
