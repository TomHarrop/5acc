setwd('/home/tom/Desktop/5accessions/samplePrep/qPCR/primerValidation/')

library(ReadqPCR)
library(NormqPCR)
library(ggplot2)
library(reshape2)

## try 2nd derivative max data
data.raw <- read.table('12.12.14 TH ACT 2HK.SDM.txt', skip = 1, header = TRUE, fill = TRUE, stringsAsFactors = FALSE, comment.char = '>', sep ='\t')
names <- paste(data.raw$Name, round(1/data.raw$Standard, 0), sep = ".")
names[grepl('^neg', names)] <- 'neg'
Detector <- rep(c('LOC_Os11g06390', 'LOC_Os03g60180', 'ACT1'), each = length(names)/3)
data.cleaned <- data.frame(Well = data.raw$Pos,PlateID = 'ACT2HK', Sample = names, Detector = Detector, Cq = data.raw$Cp, stringsAsFactors = FALSE)

data.raw.2 <- read.table('16.12.14 TH C2_LOC LHS1 APO2.txt', skip = 1, header = TRUE, fill = TRUE, stringsAsFactors = FALSE, comment.char = '>', sep = '\t')[,1:7]
data.raw.2 <- data.raw.2[data.raw.2$Include == 'True',]
names.2 <- paste(data.raw.2$Name, 1/data.raw.2$Standard, sep = ".")
names.2[grepl('^neg', names.2)] <- 'neg'
Detector.2 <- rep(c('LOC_Os08g31980', 'LHS1', 'APO2'), each = length(names.2)/3)
data.cleaned.2 <- data.frame(Well = data.raw.2$Pos,PlateID = 'C2_LOCLHS1APO2', Sample = names.2, Detector = Detector.2, Cq = data.raw.2$Cp, stringsAsFactors = FALSE)

write.table(rbind(data.cleaned, data.cleaned.2), sep = '\t', file = 'combined.cleaned.txt', quote = FALSE, row.names = FALSE)

#ACT2HK <- read.qPCR('12.12.14 TH ACT 2HK.Cp.cleaned.txt')
combined <- read.qPCR('combined.cleaned.txt')
exprs(combined)[grepl('^LHS1', rownames(exprs(combined))),'NIP.1'][c(1,2)] <- NA
ACT2HK <- combined

exprs(ACT2HK)[exprs(ACT2HK) < 18 | exprs(ACT2HK) > 36] <- NA

ACT2HK.combined <- combineTechReps(ACT2HK)
data <- exprs(ACT2HK.combined)
#data <- rbind(data, dLOC_Os03g60180 = data['LOC_Os03g60180',] - data['ACT1',], dLOC_Os11g06390 = data['LOC_Os11g06390',] - data['ACT1',])

data <- rbind(data, dAPO2 = data['APO2',] - data['ACT1',], dLHS1 = data['LHS1',] - data['ACT1',], dLOC_Os03g60180 = data['LOC_Os03g60180',] - data['ACT1',], 'dLOC_Os08g31980' = data['LOC_Os08g31980',] - data['ACT1',], 'dLOC_Os11g06390' = data['LOC_Os11g06390',] - data['ACT1',])

plotData <- melt(data, varnames = c('Detector', 'Sample'), value.name = 'Cq', na.rm = TRUE)
plotData$Species <- factor(sapply(plotData$Sample, function(x) strsplit(as.character(x), split = '\\.')[[1]][1]), levels = c('NIP', 'IR64', 'OR', 'OG', 'OB'))
plotData$Log4Concentration <- sapply(plotData$Sample, function(x) log(as.numeric(strsplit(as.character(x), split = '\\.')[[1]][2]), base = 4))


plotData2.tmp <- subset(plotData, grepl('^d', plotData$Detector))

plotData2 <- droplevels(plotData2.tmp[!plotData2.tmp$Log4Concentration > 3,])
plotData2 <- droplevels(plotData2[!(plotData2$Sample=='neg' | is.na(plotData2$Sample)),])
plotData2$Detector <- factor(plotData2$Detector, levels = c('dLOC_Os03g60180', 'dLOC_Os11g06390', 'dLHS1', 'dAPO2', 'dLOC_Os08g31980'))
#plotData2 <- plotData2[order(plotData2$Detector, plotData2$Species),]

annText <- unique(plotData2[,c('Species','Detector')])
#annText$label <- as.character(round(sapply(regressions, function(x) coef(x)[2]), 2))
#i <- 2
for (i in 1:nrow(annText)){
  s <- as.character(annText$Species[i])
  d <- as.character(annText$Detector[i])
  subData <- plotData2[plotData2$Species == s & plotData2$Detector == d,]
  label <- coef(lm(Cq ~ Log4Concentration, data = subData))['Log4Concentration']
  annText$label[i] <- paste("grad = ", round(label, 2), sep = '')
  subData.2 <- data[substring(d,2),]
  annText$labMin[i] <- round(subData.2[paste(s, '.1', sep = '')], 1)
  annText$labMax[i] <- round(subData.2[paste(s, '.64', sep = '')], 1)
  annText$lPos[i] <- subset(plotData2, Detector == d & Species == s  & Log4Concentration == 0)$Cq
  annText$rPos[i] <- subset(plotData2, Detector == d & Species == s  & Log4Concentration == 3)$Cq
  annText$midPos[i] <- mean(annText$lPos[i], annText$rPos[i])
}

annText[annText$Detector == 'dLOC_Os08g31980',c('lPos', 'rPos', 'midPos')] <- annText[annText$Detector == 'dLOC_Os08g31980',c('lPos', 'rPos', 'midPos')] + 5
annText[annText$Detector == 'dLOC_Os11g06390',c('lPos', 'rPos', 'midPos')] <- annText[annText$Detector == 'dLOC_Os11g06390',c('lPos', 'rPos', 'midPos')] + 5


g <- ggplot(data = plotData2, aes(y = Cq, x = Log4Concentration))
g + theme_minimal(base_size = 10) + 
  geom_point() +
  facet_grid(Detector ~ Species) +
  stat_smooth(method = lm, level = FALSE) +
  geom_text(data = annText, aes(x = 1.5, y = midPos, label = label), vjust = 2, size = 4) +
  geom_text(data = annText, aes(x = 0, y = lPos, label = labMin), vjust = 2, hjust = 0, size = 4) +
  geom_text(data = annText, aes(x = 3, y = rPos, label = labMax), vjust = 2, hjust = 1, size = 4)

ggsave('overview.pdf', height = 170, width = 257, units = 'mm')
