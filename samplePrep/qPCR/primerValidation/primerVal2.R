setwd('/home/tom/Desktop/5accessions/samplePrep/qPCR/primerValidation/')

library(ReadqPCR)
library(NormqPCR)
library(ggplot2)
library(reshape2)

data.raw <- read.table('12.12.14 TH ACT 2HK.Cp.txt', skip = 1, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
suffix <- rep(rep(c('0','1','2','3','4'), each = 3), 15)
names <- data.raw$Name
names[!names == 'neg'] <- paste(names[!names == 'neg'], suffix, sep = '.')
Detector <- rep(c('LOC_Os11g06390', 'LOC_Os03g60180', 'ACT1'), each = length(names)/3)
data.cleaned <- data.frame(Well = data.raw$Pos,PlateID = 'ACT2HK', Sample = names, Detector = Detector, Cq = data.raw$Cp, stringsAsFactors = FALSE)
write.table(data.cleaned, sep = '\t', file = '12.12.14 TH ACT 2HK.Cp.cleaned.txt', quote = FALSE, row.names = FALSE)

ACT2HK <- read.qPCR('12.12.14 TH ACT 2HK.Cp.cleaned.txt')

exprs(ACT2HK)[exprs(ACT2HK) < 18 | exprs(ACT2HK) > 36] <- NA

ACT2HK.combined <- combineTechReps(ACT2HK)
data <- exprs(ACT2HK.combined)
#data <- rbind(data, dLOC_Os03g60180 = data['LOC_Os03g60180',] - data['LOC_Os11g06390',])
data <- rbind(data, dLOC_Os03g60180 = data['LOC_Os03g60180',] - data['ACT1',], dLOC_Os11g06390 = data['LOC_Os11g06390',] - data['ACT1',])

plotData <- melt(data, varnames = c('Detector', 'Sample'), value.name = 'Cq', na.rm = TRUE)
plotData$Species <- factor(sapply(plotData$Sample, function(x) strsplit(as.character(x), split = '\\.')[[1]][1]), levels = c('NIP', 'IR64', 'OR', 'OG', 'OB'))
plotData$Log4Concentration <- sapply(plotData$Sample, function(x) as.numeric(strsplit(as.character(x), split = '\\.')[[1]][2]))

plotData2.tmp <- subset(plotData, grepl('^d', plotData$Detector))

plotData2 <- droplevels(plotData2.tmp[!plotData2.tmp$Log4Concentration > 3,])

regressions <- by(data = plotData2, INDICES = list(plotData2$Species, plotData2$Detector), function(x) lm(Cq ~ Log4Concentration, x))

annText <- unique(plotData2[,c('Species','Detector')])
#annText$label <- as.character(round(sapply(regressions, function(x) coef(x)[2]), 2))
#i <- 2
annText$label <- ''
for (i in 1:nrow(annText)){
  s <- as.character(annText$Species[i])
  d <- as.character(annText$Detector[i])
  data <- plotData2[plotData2$Species == s & plotData2$Detector == d,]
  label <- coef(lm(Cq ~ Log4Concentration, data = data))['Log4Concentration']
  annText$label[i] <- round(label, 2)
}


g <- ggplot(data = plotData2, aes(y = Cq, x = Log4Concentration))
g + theme_minimal() + 
  geom_point() +
  facet_wrap(Detector ~ Species, ncol = 3) +
  stat_smooth(method = lm, level = FALSE) +
  geom_text(data = annText, aes(x = 2, y = 4, label = label))

ggsave('ACT2HK.pdf', width = 14.8, height = 21.0, units = 'cm')
