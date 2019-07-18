#!/usr/bin/env Rscript

checkExtension_Synapse <- function() {
  
  # # Bug testing
  # synapseID <- "syn11581335"
  
  mapFile <- "Code/synapse.csv"
  map <- read.csv(mapFile, stringsAsFactors = F)
  colnames(map) <- c('filename','synID')
  
  
  
  # # file <- map$filename[which(map$synID == synapseID)]
  # if (! length(file)) {
  #   stop("ERROR: Unable to find synapse ID within mapping document ... recheck getSynapse.py for bugs!")
  # }
  
  files <- gsub(map$filename, pattern = ".gz", replacement = "")
  fileType <- rep(NA, length(files))
  for (i in 1:length(files)) {
    file <- files[i]
    fastq <- grep(x = file, pattern = ".fastq$", value = T)
    fq <- grep(x = file, pattern = ".fq$", value = T)
    bam <- grep(x = file, pattern = ".bam$", value = T)
    res <- ifelse(length(bam), no = ifelse(length(fastq), 
                                           no = ifelse(length(fq),
                                                       yes = "fastq", 
                                                       no = stop("File ending not recognized")),
                                           yes = "fastq"), 
                  yes = "bam")
    fileType[i] <- res
    
  }
  
  
  if (all(fileType == unique(fileType)) & unique(fileType) %in% c("fastq", "fq")) {
    # Check if you can mark mate pairs... 
    map$mate <- 1
    if (length(grep(map$filename, pattern = ".+2\\.f.+$" ))) {
      good <- length(grep(map$filename, pattern = ".+2\\.f.+$" )) == length(grep(map$filename, pattern = ".+1\\.f.+$" ))
      if (! good) {
        stop("ERROR: CANNOT DETERMINE MATE IDENTITY MANUALLY")
      }
      map$mate[grep(map$filename, pattern = ".+2\\.f.+$" )] <- 2
    } else {
      #warning("Detecting ONLY single-end reads. If this is an error, set mate identity manually... ")
      singleReads <- T
    }
    map$gz <- F
    map$gz[grep(map$filename, pattern = ".gz$")] <- T
    # Save with updated mate info
    write.csv(map, file = "Code/synapse.csv", quote = F, row.names = F)
  } else if (! all(fileType == unique(fileType))) {
    stop("ERROR: list contains a mix of file types. This is not currently supported.")
  }
  res <- unique(fileType)
  return(res)
}


# Return result
result <- checkExtension_Synapse()
cat(result)
