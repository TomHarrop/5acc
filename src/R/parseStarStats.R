#!/usr/bin/Rscript

library(data.table)

# find the most recent cutadapt output
outputDirs <- list.dirs(path = 'output', full.names = TRUE, recursive = FALSE)
cutadaptDir <- rev(sort(outputDirs[grep('cutadapt', outputDirs)]))[1]

# find the most recent STAR output
outputDirs <- dir(path = cutadaptDir, pattern = "STAR", full.names = TRUE)
starDir <- rev(sort(outputDirs[grep('STAR', outputDirs)]))[1]

# get all the final.out files
starLogFiles <- dir(path = starDir, pattern = "Log.final.out", full.names = TRUE)
starLogs <- lapply(starLogFiles, read.delim, header = FALSE, sep = "|", fill = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
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

# save for excel
wb <- xlsx::createWorkbook()
sheet <- xlsx::createSheet(wb, sheetName = 'STARlogs')
xlsx::addDataFrame(starLogs.final, sheet, showNA = FALSE, row.names = FALSE)
saveWorkbook(wb, paste0("xlsx/STARLogs.combined-", Sys.Date(), ".xlsx"))



