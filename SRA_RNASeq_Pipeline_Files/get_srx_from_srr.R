#!/usr/bin/env Rscript

get_srx_from_srr <- function(SRA_file, project_dir) {
  # SRA_file <- "example_input.txt"
  # SRA_file <- "example_input_srx.txt"
  sras <- read.table(SRA_file, stringsAsFactors = F)
  sras <- sras$V1
  con <- DBI::dbConnect(RMySQL::MySQL(), user = "public-rds-user", port = 3306, dbname="bishoplabdb",
                        password='public-user-password', host="bishoplabdb.cyss3bq5juml.us-west-2.rds.amazonaws.com")
  
  if (substr(sras[1], 1, 3) == "SRR") {
    sql <- paste0("SELECT * FROM srr_srx_map WHERE Run IN ('", paste(sras, collapse = "', '"), "');")
  } else if (substr(sras[1], 1, 3) == "SRX") {
    sql <- paste0("SELECT * FROM srr_srx_map WHERE Experiment IN ('", paste(sras, collapse = "', '"), "');")
  }
  
  
  res <- DBI::dbSendQuery(con, sql)
  resdf <- DBI::dbFetch(res)
  write.table(resdf, file = file.path(project_dir, "Code/runInfo_Table.txt"), 
              quote = F, row.names = F, sep = "\t")
}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
arg <- args[1]
arg2 <- args[2]

# Return result
get_srx_from_srr(arg, arg2)


