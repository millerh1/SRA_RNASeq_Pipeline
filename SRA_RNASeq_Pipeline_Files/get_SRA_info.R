#!/usr/bin/env Rscript

get_SRA_info <- function(SRP, project_dir) {
  
  # # Bug testing
  # SRP <- "SRP045672"
  # project_dir <- "~/RNASeq_SRA_Pipeline_Project/"
  
  
  resp <- tryCatch(
    expr = {
      httr::GET(url = "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi", 
          query = list("save"="efetch","db"="sra","rettype"="runinfo", "term"=SRP))
    }, 
    error = function(cond) {
      return(NULL)
    }
  )
  
  if (! is.null(resp)) {
    if (resp$status_code != 200) {
      msg <- paste0("returned status code of ", as.character(resp$status_code), " when contacting SRA API.")
      
      stop(msg)
    }
    ct <- suppressMessages(read.csv(text = httr::content(resp, as = "text"), stringsAsFactors = F))
    ct <- ct[which(ct$LibrarySource == "TRANSCRIPTOMIC"),]
    if (! length(ct$Run)) {
      return("err")
    }
    write.table(ct, file = file.path(project_dir, "Code/runInfo_Table.txt"), 
                quote = F, row.names = F, col.names = T, sep = "\t")
    return(1)
  } else {
    return(0)
  }
  

  
}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]
arg2 <- args[2]

# Return result
res <- get_SRA_info(arg, arg2)
cat(res)

