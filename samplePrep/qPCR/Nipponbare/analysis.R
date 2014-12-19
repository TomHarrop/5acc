setwd('/home/tom/Desktop/5accessions/samplePrep/qPCR/Nipponbare')

library(ReadqPCR)
library(NormqPCR)
library(ggplot2)

############
### MUNG ###
############

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

###############
### ANALYSE ###
###############

data <- read.qPCR('LC480.munged.txt')

## manually set values above maximum validated for each cutoff to NA
exprs(data)[grepl('LOC_Os11g06390',rownames(exprs(data))),][exprs(data)[grepl('LOC_Os11g06390',rownames(exprs(data))),] > 27.7] <- NA
exprs(data)[grepl('LOC_Os03g60180',rownames(exprs(data))),][exprs(data)[grepl('LOC_Os03g60180',rownames(exprs(data))),] > 32.5] <- NA
exprs(data)[grepl('ACT1',rownames(exprs(data))),][exprs(data)[grepl('ACT1',rownames(exprs(data))),] > 35] <- NA
exprs(data)[grepl('LHS1',rownames(exprs(data))),][exprs(data)[grepl('LHS1',rownames(exprs(data))),] > 35] <- Inf

# make case-control matrix to deal with borderline calls
stageMatrix <- cbind(c(rep(1,3), rep(0,9)), c(rep(1,3), rep(0,9))[c(10:12,1:9)], c(rep(1,3), rep(0,9))[c(7:12,1:6)], c(rep(1,3), rep(0,9))[c(4:12,1:3)])
colnames(stageMatrix) <- c("N1", "N2", "N3", "N4")
rownames(stageMatrix) <- sampleNames(data)
sMaxM <- t(as.matrix(c(2,2,2,2)))
colnames(sMaxM) <- c("N1", "N2", "N3", "N4")

# deal with the borderline Ct values
data <- makeAllNewVal(qPCRBatch = data, contrastM = stageMatrix, sampleMaxM = sMaxM, newVal = NA)

# combine technical replicates
data <- combineTechReps(data)

# select housekeepers. wtf?
res <- selectHKs(data, method = 'geNorm', minNrHKs = 3, trace = FALSE, log = FALSE, Symbols = featureNames(data))
ranks <- data.frame(c(1, 1:3), res$ranking)
names(ranks) <- c('rank', 'gene')

# delta-delta-Ct with three housekeepers
hkgs <- c('LOC_Os11g06390', 'LOC_Os03g60180', 'ACT1')
results.N1 <- deltaDeltaCq(qPCRBatch = data, hkgs = hkgs, contrastM = stageMatrix, case = "N1", control = "N4", hkgCalc = 'geom')
results.N2 <- deltaDeltaCq(qPCRBatch = data, hkgs = hkgs, contrastM = stageMatrix, case = "N2", control = "N4", hkgCalc = 'geom')
results.N3 <- deltaDeltaCq(qPCRBatch = data, hkgs = hkgs, contrastM = stageMatrix, case = "N3", control = "N4", hkgCalc = 'geom')
results.N4 <- deltaDeltaCq(qPCRBatch = data, hkgs = hkgs, contrastM = stageMatrix, case = "N4", control = "N4", hkgCalc = 'geom')

# plot
resultsList <- lapply(list(results.N1,results.N2,results.N3, results.N4), function(x) x[2,c(6:8)])
names(resultsList) <- c("N1", 'N2', 'N3', 'N4')
plotData <- do.call(rbind, resultsList)
plotData <- data.frame(apply(plotData, 2, function(x) as.numeric(x)))
names(plotData) <- c('value', 'min', 'max')
plotData$Stage <- factor(c("Rachis Meristem", 'Branch Meristem', 'Spikelet Meristem', 'Floret Meristem'), levels = c("Rachis Meristem", 'Branch Meristem', 'Spikelet Meristem', 'Floret Meristem'))

g <- ggplot(plotData, aes(x = Stage, y = value, ymin = min, ymax = max))
g + theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_bar(stat = 'identity', fill = 'grey', colour = 'black') + geom_errorbar(width = 0.2)+
  xlab(NULL) + ylab(expression(atop('Relative expression normalised to 3 HKs',(2^{- Delta * Delta * C[t]} %+-% SD))))

ggsave('LHS1_Nipponbare.pdf', height = 170*2/3, width = 257*2/3, units = 'mm')
