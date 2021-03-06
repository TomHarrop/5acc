setwd('/home/tom/Desktop/5accessions/samplePrep/qPCR/Nipponbare')

library(ReadqPCR)
library(NormqPCR)
library(ggplot2)
library(scales)

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
exprs(data)[grepl('LOC_Os11g06390',rownames(exprs(data))),][exprs(data)[grepl('LOC_Os11g06390',rownames(exprs(data))),] > 27.7] <- Inf
exprs(data)[grepl('LOC_Os03g60180',rownames(exprs(data))),][exprs(data)[grepl('LOC_Os03g60180',rownames(exprs(data))),] > 32.5] <- Inf
exprs(data)[grepl('ACT1',rownames(exprs(data))),][exprs(data)[grepl('ACT1',rownames(exprs(data))),] > 35] <- Inf
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

# calculate CIs
getCIs <- function(x){
  resSub <- subset(x, ID == 'LHS1', select = c('2^-ddCt', '2^-ddCt.min', '2^-ddCt.max'))
  resSub <- apply(resSub, 2, as.numeric)
  resSub[is.na(resSub)] <- 0
  CImax <- (1.96 * resSub["2^-ddCt.max"]) - (0.96 * resSub["2^-ddCt"])
  CImin <- (1.96 * resSub["2^-ddCt.min"]) - (0.96 * resSub["2^-ddCt"])
  output <- c(resSub, CImin=CImin, CImax=CImax)
  names(output) <- c('2^-ddCt', '2^-ddCt.min', '2^-ddCt.max', 'CI95min', 'CI95max')
  return(output)
}

resultsList <- lapply(list(results.N1,results.N2,results.N3, results.N4), getCIs)
names(resultsList) <- c("N1", 'N2', 'N3', 'N4')

# t-tests
t.testFromMeans <- function(x1, x2, sd1, sd2, n1, n2) {
  df <- n1 + n2 -2
  poolvar <- (((n1-1) * sd1^2) + ((n2 - 1) * sd2^2))/df
  t <- (x1 - x2) / sqrt(poolvar * ((1 / n1) + (1 / n2)))
  sig <- 2 * (1 - (pt(abs(t), df)))
  output <- c(t, sig, df)
  names(output) <- c('t', 'p', 'df')
  return(output)
}

plotData <- do.call(rbind, resultsList)
plotData <- data.frame(apply(plotData, 2, function(x) as.numeric(x)))
names(plotData) <- c('value', 'min', 'max', 'CI95min', 'CI95max')
plotData$Stage <- factor(c("Rachis\nMeristem", 'Branch\nMeristem', 'Spikelet\nMeristem', 'Floret\nMeristem'), levels = c("Rachis\nMeristem", 'Branch\nMeristem', 'Spikelet\nMeristem', 'Floret\nMeristem'))

plotData$p <- c(0, sapply(c(2:4), function(i) t.testFromMeans(
  resultsList[[i]]["2^-ddCt"], # x1
  resultsList[[i - 1]]["2^-ddCt"], # x2
  resultsList[[i]]["2^-ddCt"] - resultsList[[i]]["2^-ddCt.min"], # sd1
  resultsList[[i - 1]]["2^-ddCt.max"] - resultsList[[i - 1]]["2^-ddCt"], # sd2
  3,3)['p']))
plotData$eqn <- paste("italic(p) == ", gsub("e", "%*% 10 ^", scientific(plotData$p, digits = 3)), sep = ' ')

plotData['eqn'][!plotData['p'] > 0] <- NA
plotData$eqn[is.nan(plotData$p)] <- NA

g <- ggplot(plotData, aes(x = Stage, y = value, ymin = CI95min, ymax = CI95max))
g <- g + theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_bar(stat = 'identity', fill = 'grey', colour = 'black') + geom_errorbar(width = 0.2) +
  xlab(NULL) + ylab(expression(atop('Normalised relative expression',(2^{- Delta * Delta * C[t]} %+-% "95% CI, " *italic(n) == 3 ))))

# add the annotated p-value. Could 'dodge' it to the left to get p-value label on line segment
g <- g + geom_text(aes(y = CI95max+0.05, label = eqn), size = 3, parse = TRUE, angle = 90, vjust = 0.5, hjust = 0) + ylim(c(0,2))



g + annotate('segment', x = 2, xend = 3, y = plotData[3,'CI95max'] + 0.1, yend = plotData[3,'CI95max'] + 0.1, colour = 'black') +
  annotate('segment', x = 3, xend = 3, y = plotData[3,'CI95max'] + 0.05, yend = plotData[3,'CI95max'] + 0.1, colour = 'black') +
  annotate('segment', x = 2, xend = 2, y = plotData[2,'CI95max'] + 0.05, yend = plotData[3,'CI95max'] + 0.1, colour = 'black')

plotData[3, 'p']

ggsave('LHS1_Nipponbare.pdf', height = 170*2/3, width = 257*2/4, units = 'mm')
