setwd('/home/tom/Desktop/5accessions/samplePrep/qPCR/Nipponbare')

library(ReadqPCR)
library(NormqPCR)

LC480.absQuant <- read.table('TH 18.12.14 Nipponbare RQs.txt', skip = 1, header = TRUE, fill = TRUE, stringsAsFactors = FALSE, comment.char = '>', sep = '\t')[,1:7]

LOC_Os11g06390.positions <- c(paste('A', seq(1,18), sep = ''), paste('E', seq(1,18), sep = ''))
LOC_Os03g60180.positions <- c(paste('B', seq(1,18), sep = ''), paste('F', seq(1,18), sep = ''))
ACT1.positions <- c(paste('C', seq(1,18), sep = ''), paste('G', seq(1,18), sep = ''))
LHS1.positions <- c(paste('D', seq(1,18), sep = ''), paste('H', seq(1,18), sep = ''))

LC480.absQuant$Detector <- ''
LC480.absQuant[LC480.absQuant$Pos %in% LOC_Os11g06390.positions,]$Detector <- 'LOC_Os11g06390'
LC480.absQuant[LC480.absQuant$Pos %in% LOC_Os03g60180.positions,]$Detector <- 'LOC_Os03g60180'
LC480.absQuant[LC480.absQuant$Pos %in% ACT1.positions,]$Detector <- 'ACT1'
LC480.absQuant[LC480.absQuant$Pos %in% LHS1.positions,]$Detector <- 'LHS1'

LC480.absQuant.munged <- data.frame(Well = LC480.absQuant$Pos,PlateID = 'get_plate_ID', Sample = LC480.absQuant$Name, Detector = LC480.absQuant$Detector, Cq = LC480.absQuant$Cp, stringsAsFactors = FALSE)
LC480.absQuant.munged <- LC480.absQuant.munged[!(grepl('neg', LC480.absQuant.munged$Sample, ignore.case = TRUE) | grepl('sample', LC480.absQuant.munged$Sample, ignore.case = TRUE) | grepl('cal', LC480.absQuant.munged$Sample, ignore.case = TRUE)),]
write.table(LC480.absQuant.munged, sep = '\t', file = 'LC480.munged.txt', quote = FALSE, row.names = FALSE)

data <- read.qPCR('LC480.munged.txt')
exprs(data)[grepl('LHS1',rownames(exprs(data))),]

data2 <- combineTechReps(data)

stageMatrix <- cbind(c(rep(1,3), rep(0,9)), c(rep(1,3), rep(0,9))[c(10:12,1:9)], c(rep(1,3), rep(0,9))[c(7:12,1:6)], c(rep(1,3), rep(0,9))[c(4:12,1:3)])
colnames(stageMatrix) <- c("N1", "N2", "N3", "N4")
rownames(stageMatrix) <- sampleNames(data2)
hkgs <- c('LOC_Os11g06390', 'LOC_Os03g60180')
results <- deltaDeltaCq(qPCRBatch = data, hkgs = 'ACT1', contrastM = stageMatrix)
