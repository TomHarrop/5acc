#!/usr/bin/env Rscript

library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}
GenerateMessage(
  "Download GWAS data from Crowell et al. 10.1038/ncomms10527.")

# download the data from crowell paper
GenerateMessage("Requesting tar file from ricediversity.org")
crowell.url <- "https://ricediversity.org/data/sams_data/GWAS-And-Association-Network-Results-Suppl-1-11.tar.gz"
temp <- tempfile(fileext = ".tar.gz")
download.file(crowell.url, temp)

# where to untar files
tar.dir <- tempdir()

# binned GWAS results
binned.gwas.file <- "GWAS-And-Association-Network-Results-Suppl-1-11/Supplementary_Data_3/GWAS_Results_NetworkData_HD_Covariate.txt"

GenerateMessage("Unpacking Supplementary_Data_3")
untar(temp, files = paste0("./", binned.gwas.file), exdir = tar.dir)

GenerateMessage("Generating R data")
binned.gwas.crowell <- data.table(
  read.table(
    paste(tar.dir, binned.gwas.file, sep = "/"),
    sep = "\t", header = TRUE, stringsAsFactors = FALSE))

# save output
out.dir <- "data/gwas/crowell/"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(binned.gwas.crowell, paste0(out.dir, "/binned_gwas_crowell.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)
