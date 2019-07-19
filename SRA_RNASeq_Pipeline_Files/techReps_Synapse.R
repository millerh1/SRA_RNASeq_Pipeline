#!/usr/bin/env Rscript

techReps_Synapse <- function(synapseID, mergeFile) {
  
  # # # # # Bug testing
  # mergeFile <- "~/bin/SRA_RNASeq_Pipeline_Files/synapseMergeFileExample.txt"
  # synapseID <- "syn11581335"
  # mapFile <- "~/Bishop.lab/Preprocessing/RNA_Seq/ROSMAP_microGlia/Code/synapse.csv"
  
  #cat("\nGetting merge info for", synapseID, "...\n")
  #cat(' ')
  mapFile <- "Code/synapse.csv"
  map <- read.csv(mapFile, stringsAsFactors = F)
  mergeFile <- read.table(mergeFile, sep = "\t", fill = T, header = T, stringsAsFactors = F)
  colnames(mergeFile) <- c("filename", "sampleID")
  colnames(map) <- c('filename','synID', "mate", "gz")
  mergeFile <- merge(x = map, y = mergeFile, by = "filename", all.x = T)
  mergeFile <- mergeFile[order(mergeFile$sampleID, mergeFile$mate, mergeFile$synID),]
  # Now the mergeFile is in order for concatenation
  mergeFile$sampleID[which(is.null(mergeFile$sampleID))] <- mergeFile$synID[which(is.null(mergeFile$sampleID))]
  mergeFile$newName <- paste0(mergeFile$sampleID, "_", mergeFile$mate, ".fastq")

  if (synapseID %in% mergeFile$synID) {
    sampleID <- mergeFile$sampleID[which(mergeFile$synID == synapseID)]
    # Check to make sure files don't already exist
    newNames <- mergeFile$newName[which(mergeFile$sampleID == sampleID)]
    newNames <- unique(newNames)
    if (all(file.exists(file.path("Data/Raw_Reads", newNames)))) {
      res <- sampleID
      return(res)
    }
    
    files <- mergeFile$synID[which(mergeFile$sampleID == sampleID)]
    ind <- which(files == synapseID)
    
    if (ind == 1) {
      sampleInfo <- mergeFile[which(mergeFile$sampleID == sampleID),]
      sampleInfo1 <- sampleInfo[which(sampleInfo$mate == 1),]
      for (i in 1:length(sampleInfo1$synID)) {
        synID <- sampleInfo1$synID[i]
        # Get the data from synapse

        cmd <- paste0("synapse get ", synID," --downloadLocation Data/Raw_Reads/synDownloadTmp/")
        system(cmd, ignore.stdout = T)
        # Gunzip if applicable
        gz <- sampleInfo1$gz[i]
        downFile <- sampleInfo1$filename[i]
        if (gz) {
          cmd <- paste0("gunzip Data/Raw_Reads/synDownloadTmp/", downFile)
          system(cmd, ignore.stdout = T)
          downFile <- gsub(downFile, pattern = ".gz",replacement =  "")
        }
        # Concattenate the files
        newName <- sampleInfo1$newName[i]
        cmd <- paste0("cat Data/Raw_Reads/synDownloadTmp/", downFile,
                      " >> ", file.path("Data/Raw_Reads", newName))
        system(cmd)
        # Delete raw files
        system(paste0("rm -rf Data/Raw_Reads/synDownloadTmp/", downFile), ignore.stdout = T)

      }
      sampleInfo2 <- sampleInfo[which(sampleInfo$mate == 2),]
      if (length(sampleInfo2)) {
        for (i in 1:length(sampleInfo2$synID)) {
          synID <- sampleInfo2$synID[i]
          # Get the data from synapse
          newName <- sampleInfo2$newName[i]
          cmd <- paste0("synapse get ", synID," --downloadLocation Data/Raw_Reads/synDownloadTmp/")
          system(cmd, ignore.stdout = T)
          # Gunzip if applicable
          gz <- sampleInfo2$gz[i]
          downFile <- sampleInfo2$filename[i]
          if (gz) {
            cmd <- paste0("gunzip Data/Raw_Reads/synDownloadTmp/", downFile)
            system(cmd, ignore.stdout = T)
            downFile <- gsub(downFile, pattern = ".gz",replacement =  "")
          }
          # Concattenate the files
          cmd <- paste0("cat Data/Raw_Reads/synDownloadTmp/", downFile,
                        " >> ", file.path("Data/Raw_Reads", newName))
          system(cmd)

          # Delete raw files
          system(paste0("rm -rf Data/Raw_Reads/synDownloadTmp/", downFile), ignore.stdout = T)

        }
      }

      res <- sampleID
      #print(res)
      #print("RES")
    } else {
      res <- "skip"
    }
  } else {
    #msg <- paste0("In techReps_Synapse.R: could not locate synapse ID ", synapseID,
                  #" in mergeFile provided, even though it belongs within the selected synapse directory. Add it to the exlusion list or the mergeFile to run it. Skipping this accession ... ")
    #warning(msg)
    res <- "no"
  }
  
  return(res)

}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]
argTable <- args[2]


# Return result
result <- suppressMessages(techReps_Synapse(arg, argTable))
cat(result)
