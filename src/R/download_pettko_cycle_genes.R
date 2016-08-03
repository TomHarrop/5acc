#!/usr/bin/env Rscript

library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}
GenerateMessage("Retrieve cell cycle genes from 10.1007/s10265-015-0754-3")

# download rice cell cycle genes from Pettkó-Szandtner et al. 2015
# (10.1007/s10265-015-0754-3)
GenerateMessage("Downloading file")
tmp <- tempfile(fileext = ".xls")
pettko.url <- "https://static-content.springer.com/esm/art%3A10.1007%2Fs10265-015-0754-3/MediaObjects/10265_2015_754_MOESM1_ESM.xls"
download.file(pettko.url, tmp)

GenerateMessage("Tidying data")

# read in from xls
pettko.raw.df <- xlsx::read.xlsx(tmp, sheetIndex = 1, stringsAsFactors = FALSE)
pettko.raw <- data.table(pettko.raw.df)
pettko.cols <- pettko.raw[!is.na(proposed.gene.name), .(
  proposed.gene.name, rice.gene.locus.identifier
)]

# authors confused, fix identifiers
pettko.cols[, gene.id := gsub("[^[:alnum:]]", "", rice.gene.locus.identifier)]
pettko.cols[!grep("^Os", gene.id), gene.id := NA]
pettko.cols[!is.na(gene.id), gene.id := paste0("LOC_", gene.id)]

# extract gene families from proposed names
pettko.cols[!is.na(gene.id), class := gsub("^[^_]+_", "", proposed.gene.name)]
pettko.cols[!is.na(gene.id), class := gsub("_[^_]+$", "", class)]
pettko.cols[!is.na(gene.id), class := gsub("[^[:alpha:]]+$", "", class)]

# select results
setkey(pettko.cols, gene.id)
pettko.genes <- unique(pettko.cols[!is.na(gene.id), .(gene.id, class)])

# save output
out.dir <- "data/goi"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(pettko.genes, paste0(out.dir, "/pettko_cycle.genes.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.pettko_cycle.genes.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)
