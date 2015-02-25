#!/usr/bin/Rscript

library(gdata)
library(reshape2)
library(xlsx)

setwd('/home/tom/Desktop/5accessions/samplePrep/nanoQuant/')

sheets <- c('OR', 'OG', 'OB', 'IR64')
colClasses <- c('NULL', rep('character', 10), 'NULL')
names <- paste(rep(paste('N', seq(1:4), sep = ''), each = 3), paste('R', seq(1:3), sep = ''), sep = '')

getTidyData <- function(sheet){
  rawData <- read.xls(xls = 'NQ20.11.14.xls', sheet = sheet, skip = 2, nrows = 23, colClasses = colClasses, header = FALSE)
  rawData <- rawData[!rawData$V3 == 'Abs',]
  toNum <- colnames(rawData)[c(1,2,4,6,7,9)]
  df.tmp <- data.frame(nrows = 16) 
  for (col in toNum){
    newCol <- as.numeric(rawData[,col])
    df.tmp <- cbind(df.tmp, newCol)
  }
  df.tmp <- subset(df.tmp, select = -c(nrows))
  rawData.split1 <- df.tmp[,c(1:3)]
  rawData.split2 <- df.tmp[c(1:8),c(4:6)]
  names(rawData.split2) <- names(rawData.split1)
  tidyData <- rbind(rawData.split1, rawData.split2)
  colnames(tidyData) <- c('wl', 'reading', 'value')
  tidyData$sample <- rep(names, each = 2)
  tidyData$type <- factor(rep(c('concentration', 'ratio'), times = 12))
  dcast(tidyData, sample ~ type)
}

nanoQuantData.list <- lapply(sheets, getTidyData)
names(nanoQuantData.list) <- sheets
nanoQuantData <- subset(do.call(cbind, nanoQuantData.list), select = -c(OG.sample, OB.sample, IR64.sample))
rownames(nanoQuantData) <- nanoQuantData$OR.sample
nanoQuantData <- subset(nanoQuantData, select = -c(OR.sample))

nanoQuantData <- nanoQuantData[,c(7,8,1:6)]

# add repeated IR64 N3R1 data
nanoQuantData['N3R1',c('IR64.concentration', 'IR64.ratio')] <- c(94.7,1.86)
# add Nipponbare data (done on nanoDrop)
nanoQuantData[,'NIP.concentration'] <- c(79, 376.8, 389.6, 357.9, 347, 309.8, 586.6, 411.5, 709, 457.1, 530.1, 510.9)
nanoQuantData[,'NIP.ratio'] <- c(2.02, 1.99, 1.99, 2.01, 1.99, 1.99, 2.01, 1.99, 2.01, 1.95, 2.02, 2.01)
nanoQuantData <- nanoQuantData[,c('NIP.concentration', 'NIP.ratio', 'IR64.concentration', 'IR64.ratio', 'OR.concentration', 'OR.ratio', 'OG.concentration', 'OG.ratio', 'OB.concentration', 'OB.ratio')]
# add repeated Nipponbare data (N1R1)
nanoQuantData['N1R1',c('NIP.concentration', 'NIP.ratio')] <- c(148.0,1.87)

saveRDS(nanoQuantData,file = 'nanoQuantData.RDS')
#nanoQuantData <- readRDS('nanoQuantData.RDS')

colnames(nanoQuantData) <- gsub('\\.', '\n', colnames(nanoQuantData))
colnames(nanoQuantData) <- gsub('\nratio', '\n260/280', colnames(nanoQuantData))

write.xlsx(nanoQuantData, file = 'nanoQuantData.xlsx', sheetName = 'Sheet1')

###########################
### LIBRARY PREP OUTPUT ###
###########################

libPreps <- nanoQuantData[grepl('^N[2|3]', rownames(nanoQuantData)), grepl('conc', colnames(nanoQuantData))]
volLib <- apply(libPreps, 2, function(x) 400/x)
colnames(volLib) <- gsub("concentration", "vol", colnames(volLib))
volH2O <- apply(volLib, 2, function(x) 10-x)
colnames(volH2O) <- gsub("vol", "H2O", colnames(volH2O))

output <- data.frame(libPreps, volLib, volH2O)
#paste(colnames(output)[rep(x, 5) + rep(c(0,1,2,3,4), each = 3)], collapse = "', '")
output <- output[,c('NIP.concentration', 'NIP.vol', 'NIP.H2O', 'IR64.concentration', 'IR64.vol', 'IR64.H2O', 'OR.concentration', 'OR.vol', 'OR.H2O', 'OG.concentration', 'OG.vol', 'OG.H2O', 'OB.concentration', 'OB.vol', 'OB.H2O')]
colnames(output) <- gsub('\\.', '\n', colnames(output))

write.xlsx(round(output, 1), file = 'libPreps.xlsx', sheetName = 'libPreps')
