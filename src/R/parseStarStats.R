#!/usr/bin/Rscript

library(data.table)

# find the Log.final.out files
starLogFiles <- list.files(path = 'output', pattern = "Log.final.out",
                           full.names = TRUE, recursive = TRUE)
starLogFiles <- starLogFiles[!grepl("step1", starLogFiles)]
starLogFiles <- starLogFiles[!grepl("remap", starLogFiles)]
starLogs <- lapply(starLogFiles, read.delim, header = FALSE, sep = "|", fill = TRUE,
                   strip.white = TRUE, stringsAsFactors = FALSE)
names(starLogs) <- gsub('.Log.final.out', '', basename(starLogFiles), fixed = TRUE)

# combine the logs into one data.frame
starLogs.frame <- do.call(cbind, args = lapply(starLogs, function(x) x[2]))
data.table::setnames(starLogs.frame, names(starLogs))
rownames(starLogs.frame) <- c(starLogs[[1]]$V1)

# transpose
starLogs.frame.t <- t(starLogs.frame)

# get the columns into the right format..
starLogs.dt <- as.data.table(starLogs.frame.t, keep.rownames = TRUE)
setnames(starLogs.dt, 'rn', 'Library')

# 1. remove empty columns
starLogs.dt[,c("UNIQUE READS:", "MULTI-MAPPING READS:", "UNMAPPED READS:") := NULL]

# 2. fix columns that are actually in percent
setnames(starLogs.dt, old = c("Deletion rate per base", "Insertion rate per base"),
         new = c("Deletion rate per base, %", "Insertion rate per base, %"))

# 3. remove percent signs and convert numeric columns
tonumeric <- names(starLogs.dt[, !c("Library", "Started job on", "Started mapping on", "Finished on"), with = FALSE])
stripPercent <- function(x){gsub("%", "", x, fixed = TRUE)}
starLogs.dt.noPercent <- starLogs.dt[,lapply(.SD, stripPercent), .SDcols = tonumeric]
starLogs.dt.num <- starLogs.dt.noPercent[,lapply(.SD, as.numeric), .SDcols = tonumeric]

# 4. convert date columns to POSIX
dateConvert <- function(x) {as.POSIXct(strptime(x, format = '%b %d %H:%M:%S'))}
todate <- names(starLogs.dt[, c("Started job on", "Started mapping on", "Finished on"), with = FALSE])
starLogs.date <- starLogs.dt[, lapply(.SD, dateConvert), .SDcols = todate]

# 5. cbind it together
starLogs.final <- cbind(starLogs.dt[,Library], starLogs.date, starLogs.dt.num)
setnames(starLogs.final, 'V1', 'Library')
setkey(starLogs.final, "Library")

# 6. get the number of reads mapped to the Nipponbare GTF
readsPerGeneFiles <- list.files(path = 'output', pattern = "ReadsPerGene.out.tab",
                           full.names = TRUE, recursive = TRUE)
readsPerGene <- lapply(readsPerGeneFiles, read.table, stringsAsFactors = FALSE)
names(readsPerGene) <- gsub(".*([[:upper:]][[:digit:]]).*", "\\1", readsPerGeneFiles)
tidyReads <- function(x) {
  x <- as.data.table(x)
  setnames(x, old = c("V1", "V4"), new = c("msuId", "reverse"))
  x[, c("V2", "V3") := NULL]
  names(x) <- sub("V4", "reverse", names(x))
  return(data.frame(x, row.names = "msuId"))
}
tidyReadsFrames <- lapply(readsPerGene, tidyReads)
tidyReadsFrame <- do.call(cbind, tidyReadsFrames)
setnames(tidyReadsFrame, names(tidyReadsFrames))
readsTable <- as.data.table(tidyReadsFrame, keep.rownames = TRUE)
setnames(readsTable, "rn", "msuId")
sdc <- names(readsTable)[!names(readsTable) == "msuId"]
readSums <- data.table(t(readsTable[!grep("^N", msuId), lapply(.SD, sum),
                                    .SDcols = sdc]), keep.rownames = TRUE)
setnames(readSums, c("rn", "V1"), c("Library", "Number of reads in genes"))
setkey(readSums, "Library")

starLogs.final <- merge(starLogs.final, readSums)

# MAKE OUTPUT FOLDER
outDir <- "output/mappingStats"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# save RDS
saveRDS(starLogs.final, paste0(outDir, "/starLogs.Rds"))

# save CSV
write.table(as.data.frame(starLogs.final), paste0(outDir, "/starLogs.tsv"),sep = "\t",
            quote = FALSE, na = "",
          row.names = FALSE, col.names = TRUE)

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/SessionInfo.txt")
writeLines(sInf, logLocation)

quit(save = "no", status = 0)
