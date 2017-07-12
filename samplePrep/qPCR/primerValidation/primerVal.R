setwd('/home/tom/Desktop/5accessions/samplePrep/qPCR/primerValidation/')

library(ReadqPCR)
library(NormqPCR)

cycData <- read.LC480('12.12.14 TH ACT 2HK.Cp.txt', skip = 1, header = TRUE, fill = TRUE)
sampleInfo <- read.LC480SampleInfo('sampleInfo.txt')

cycData1 <- merge(cycData, sampleInfo)
suffix <- rep(rep(c('1','4','16','64','256'), each = 3), 15)
names <- pData(cycData1)$'Sample name'
names[!names == 'neg'] <- paste(names[!names == 'neg'], suffix, sep = '.')
pData(cycData1)$'Sample name' <- names

combineTechReps(cycData1)
?read.LC480
