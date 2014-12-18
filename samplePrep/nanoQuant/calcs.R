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

saveRDS(nanoQuantData,file = 'nanoQuandData.RDS')
#nanoQuantData <- readRDS('nanoQuandData.RDS')

colnames(nanoQuantData) <- gsub('\\.', '\n', colnames(nanoQuantData))
colnames(nanoQuantData) <- gsub('\nratio', '\n260/280', colnames(nanoQuantData))

write.xlsx(nanoQuantData, file = 'nanoQuantData.xlsx', sheetName = 'Sheet1')
