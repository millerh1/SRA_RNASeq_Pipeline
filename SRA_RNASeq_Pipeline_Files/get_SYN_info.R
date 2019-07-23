#!/usr/bin/env Rscript

get_SYN_info <- function(synID, synList, scriptDIR) {
  
  # # # # # Bug testing
  # synID <- "syn4894912"
  # # synList <- "~/Bishop.lab/Preprocessing/RNA_Seq/MayoRNAseq_Full2/synList1.csv"
  # scriptDIR <- "~/bin/SRA_RNASeq_Pipeline_Files/"
  
  
  dbdat <- read.csv(file.path(scriptDIR, "synapse_studies_list.csv"), stringsAsFactors = F)
  
  dbdat2 <- NULL
  dbdat2 <- dbdat[which(dbdat$parentId == synID | dbdat$projectId == synID),]
  if (synList != 'none') {
    synList <- read.table(synList, header = F, stringsAsFactors = F)
    dbdat3 <- dbdat[which(dbdat$specimenID %in% synList[,1] | dbdat$name %in% synList[,1]),]
    dbdat2 <- dbdat3[which(dbdat3$parentId %in% dbdat2$parentId),]
  }
  colnames(dbdat2)[3] <- "Run"
  timevec <- Sys.time()
  timevec <- gsub(timevec, pattern =  ":", replacement =  ".")
  timevec <- gsub(timevec, pattern =  " ", replacement =  "_")
  fileStr <- paste0("Code/runInfo_Table_", timevec, ".txt")
  write.table(dbdat2, file = fileStr, 
              quote = F, row.names = F, col.names = T, sep = "\t")
  return(fileStr)
}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]
arg2 <- args[2]
arg3 <- args[3]

# Return result
res <- suppressWarnings(suppressMessages(get_SYN_info(arg, arg2, arg3)))
cat(res)

