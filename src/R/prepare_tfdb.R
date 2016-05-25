#!/usr/bin/env Rscript

library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

GenerateMessage("Downloading TFDB annotation")

# download tfdb to tempfile
tmp <- tempfile(fileext = ".tab")
call_download <- paste(c(
  'curl --data-urlencode "_submit_check=1"',
  '--data-urlencode "DOWNLOAD=list of transcription factors"',
  '--data-urlencode "SUBMIT=Download"',
  'http://plntfdb.bio.uni-potsdam.de/v3.0/export.php > ',
  tmp), collapse = " ")
system(call_download)

# read in
GenerateMessage("Reading file")
tfdb.raw <- data.table(read.table(tmp, header = T, sep = '\t',
                                 stringsAsFactors = F))

# tidy
GenerateMessage("Tidying data")
setkey(tfdb.raw, 'Species', 'Family', 'Protein.ID', 'Category')
tfdb.raw[, Genome.EST := NULL]
tfdb <- unique(tfdb.raw)

# save
out.dir <- "data/tfdb"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

# full table not only rice
saveRDS(tfdb, paste0(out.dir, "/tfdb.Rds"))

# rice only
tfdb.os <- tfdb[Species == "Oryza sativa subsp. japonica", .(
  Protein.ID, Family, Category
)]
saveRDS(tfdb.os, paste0(out.dir, "/tfdb_os.Rds"))

# logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.wald_stage_species.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)
