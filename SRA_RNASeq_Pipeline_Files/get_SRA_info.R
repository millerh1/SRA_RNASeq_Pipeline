#!/usr/bin/env Rscript

get_SRA_info <- function(SRP, SRAList) {
  
  # # # # # # Bug testing
  # SRP <- "SRR5527012"
  # # SRAList <- "~/Bishop.lab/EWS_CTR_All_Cells/Code/runInfo_Table_2019-07-22_13.59.15.accessionList.txt"
  # 
  
  toCheck <- c()
  if (SRP != 'none') {
    toCheck <- SRP
  }
  
  if (SRAList != 'none') {
    SRAList <- read.table(SRAList, header = F, stringsAsFactors = F)
    toCheck <- c(toCheck, SRAList$V1)
  }
  goodStr <- c("SRP", "SRR", "SRX")
  toCheck <- toCheck[which(substr(toCheck, 1, 3) %in% goodStr)]
  
  for (i in 1:length(toCheck)) {
    term <- toCheck[i]
    resp <- NULL
    while(is.null(resp)) {
      resp <- tryCatch(
        expr = {
          Sys.sleep(.5)
          httr::GET(url = "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi", 
                    query = list("save"="efetch","db"="sra","rettype"="runinfo", "term"=term))
        }, 
        error = function(cond) {
          Sys.sleep(.5)
          NULL
        }
      )
    }
    
    if (resp$status_code != 200) {
      msg <- paste0("ERROR: returned status code of ", 
                    as.character(resp$status_code),
                    " when contacting SRA API to obtain data for  ", term)
      return(msg)
    }
    ct <- suppressMessages(read.csv(text = httr::content(resp, as = "text"), stringsAsFactors = F))
    ct <- ct[which(ct$LibrarySource == "TRANSCRIPTOMIC"),]
    if (i == 1) {
      runInfo <- ct
    } else {
      cols <- which(colnames(runInfo) %in% colnames(ct))
      runInfo <- merge(x = runInfo, y = ct, by = cols, all = T)
    }
    
  }
  
  # Now get the group info and location info
  groupInfo <- c()
  locationInfo <- c()
  queryVec <- runInfo$BioSample
  for (i in 1:length(queryVec)) {
    term <- queryVec[i]
    resp <- NULL
    while(is.null(resp)) {
      resp <- tryCatch(
        expr = {
          Sys.sleep(.5)
          httr::GET(url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", 
                    query = list("db"="biosample", "id"=term))
        }, 
        error = function(cond) {
          NULL
        }
      )
    }
    if (resp$status_code != 200) {
      gi <- NA
      li <- NA
    } else {
      res <- XML::xmlToList(httr::content(resp, as = "text"))
      gi <- paste0(res$BioSample$Owner$Contacts$Contact$Name$First, " ",
                   res$BioSample$Owner$Contacts$Contact$Name$Last)
      li <- res$BioSample$Owner$Name
    }
    groupInfo <- c(groupInfo, gi)
    locationInfo <- c(locationInfo, li)
  }
  n <- length(colnames(runInfo))-2
  runInfo2 <- runInfo[,c(1, 7, 11, 13, 14, 15, 16, 17, 19, 20, 21, 22, 24, 25, 27, 28, 29:n)]
  runInfo2 <- as.data.frame(runInfo2)
  runInfo2$group <- groupInfo
  runInfo2$location <- locationInfo
  
  
  
  timevec <- Sys.time()
  timevec <- gsub(timevec, pattern =  ":", replacement =  ".")
  timevec <- gsub(timevec, pattern =  " ", replacement =  "_")
  timevec <- gsub(timevec, pattern =  "-", replacement =  ".")
  fileStr <- paste0("Code/runInfo_Table_", timevec, ".csv")
  #fileStr <- paste0("runInfo_Table_", timevec, ".csv")
  write.csv(runInfo2, file = fileStr, 
              quote = F, row.names = F, sep = ",")


  
  return(fileStr)
}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]
arg2 <- args[2]

# Return result
res <- suppressWarnings(suppressMessages(get_SRA_info(arg, arg2)))
cat(res)

