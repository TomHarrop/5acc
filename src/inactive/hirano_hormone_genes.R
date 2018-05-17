#!/usr/bin/env Rscript

library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}
GenerateMessage("Hormone-related genes from doi 10.1093/pcp/pcn123")

# download supplementary data from Hirano et al 2008  (10.1093/pcp/pcn123)
GenerateMessage("Downloading file")
hirano.url <-
  "http://pcp.oxfordjournals.org/content/suppl/2008/08/21/pcn123.DC1/pcp-2008-e-00334-File014.xls"
tmp <- tempfile(fileext = ".xls")
download.file(hirano.url, tmp)

# is this the nastiest xls ever? try to parse it
GenerateMessage("Parsing")
hirano.raw <- xlsx::read.xlsx(tmp, sheetIndex = 1, colIndex = c(2:7),
                              startRow = 3, stringsAsFactors = FALSE)
hirano.table <- data.table(hirano.raw)

# column names are in the second row now
GenerateMessage("Tidying")
real.col.names <- as.character(hirano.table[`NA.` == "Gene name"][1])
setnames(hirano.table, names(hirano.table), real.col.names)

# remove the columns that repeat the column names
hirano.table <- hirano.table[!`Gene name` == "Gene name"]

# remove footnotes
hirano.genes <- hirano.table[!grepl("^[a|b][[:space:]]", `Gene name`)]

# parse the hormone class. it's in the rows that have NA for every other column
hormone.classes <- hirano.genes[!is.na(`Gene name`) & is.na(`RAP-DB Name`),
                                `Gene name`]
hormone.idx <- hirano.genes[,which(`Gene name` %in% hormone.classes)]
row.start <- hormone.idx + 1
row.stop <- sapply(1:(length(row.start) - 1), function(i)
  hormone.idx[i + 1] - 1)
row.stop[length(row.start)] <- dim(hirano.genes)[1]
hirano.list <- lapply(1:length(hormone.classes), function(i)
  hirano.genes[row.start[i]:row.stop[i]])
names(hirano.list) <- hormone.classes
hirano.tidy <- rbindlist(hirano.list, idcol = "hormone.class")

# fix column names
setnames(hirano.tidy,
         c("Gene name", "Chromosomal Position", "RAP-DB Name", "probe on 44K"),
         c("gene.name", "chr.pos", "rap.id", "array.probe"))

# strip whitespace from rap ID
hirano.tidy[, rap.id := gsub("[[:space:]]", "", rap.id)]
hirano.tidy[rap.id == "notfound", rap.id := NA]

# convert rap-db ids to MSU ids
GenerateMessage("Adding MSU IDs")
msu.id.table <- hirano.tidy[!is.na(rap.id),
                            .(msu.id = oryzr::RapToMsu(rap.id)$msu.id),
                            by = rap.id]
setkey(msu.id.table, rap.id)
setkey(hirano.tidy, rap.id)
hormone.genes <- hirano.tidy[msu.id.table]

# add official symbols
GenerateMessage("Adding official gene symbols")
hormone.genes[!is.na(msu.id),
              official.symbol := oryzr::LocToGeneName(msu.id)$symbols,
              by = msu.id]
setcolorder(hormone.genes, c("msu.id", "rap.id", "official.symbol", "gene.name", "hormone.class", "cDNA",
                           "chr.pos", "array.probe", "Reference"))

# save output
out.dir <- "data/goi"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(hormone.genes, paste0(out.dir, "/hirano_hormone_genes.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.hirano_hormone_genes.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)
