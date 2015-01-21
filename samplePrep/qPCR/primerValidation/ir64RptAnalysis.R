library(ReadqPCR)
library(ggplot2)
library(reshape2)

setwd('/home/tom/Desktop/5accessions/samplePrep/qPCR/primerValidation/')

############
### MUNG ###
############

# read the LC480 text output
data.raw.2 <- read.table('/home/tom/Desktop/5accessions/samplePrep/qPCR/primerValidation/TH 18.12.14 IR64 rpt SCs.txt', skip = 1, header = TRUE, fill = TRUE, stringsAsFactors = FALSE, comment.char = '>', sep = '\t')[,1:7]

# add gene names
data.raw.2$Detector <- rep(c('LOC_Os11g06390', 'LOC_Os03g60180', 'ACT1', 'LHS1'), each = dim(data.raw.2)[1]/4)

# fix read error?
data.raw.2[data.raw.2$Name == 'error','Name'] <- 'IR64.9'
data.raw.2[data.raw.2$Name == 'error2','Name'] <- 'IR64.27'
data.raw.2[data.raw.2$Name == 'error3','Name'] <- 'IR64.9'
data.raw.2[data.raw.2$Name == 'error4','Name'] <- 'IR64.27'

# fix names
data.raw.2[grepl('^neg', data.raw.2$Name, ignore.case = TRUE), 'Name'] <- 'neg'

# write to text file
data.cleaned.2 <- data.frame(Well = data.raw.2$Pos, PlateID = 'IR64RPT', Sample = data.raw.2$Name, Detector = data.raw.2$Detector, Cq = data.raw.2$Cp, stringsAsFactors = FALSE)
write.table(data.cleaned.2, sep = '\t', file = 'IR64rpt.cleaned.txt', quote = FALSE, row.names = FALSE)

###############
### ANALYSE ###
###############

# read cleaned data with readqpcr
ir64rpt <- read.qPCR('IR64rpt.cleaned.txt')

# remove pipetting errors
exprs(ir64rpt)[exprs(ir64rpt) == Inf] <- NA

# combine technical replicates
ir64rpt.combined <- combineTechReps(ir64rpt)

# calculate delta-Ct
data <- exprs(ir64rpt.combined)
data <- rbind(data, 'dLHS1' = data['LHS1',] - data['ACT1',], 'dLOC_Os03g60180' = data['LOC_Os03g60180',] - data['ACT1',], 'dLOC_Os11g06390' = data['LOC_Os11g06390',] - data['ACT1',])

# prepare plot data frame
plotData <- melt(data, varnames = c('Detector', 'Sample'), value.name = 'Cq', na.rm = TRUE)
plotData$Species <- factor(sapply(plotData$Sample, function(x) strsplit(as.character(x), split = '\\.')[[1]][1]), levels = c('NIP', 'IR64', 'OR', 'OG', 'OB'))
plotData$Log3Concentration <- sapply(plotData$Sample,
                                     function(x) log(as.numeric(strsplit(as.character(x), split = '\\.')[[1]][2]), base = 3))

# subset plot data to only include deltaCt values and remove negative and
# out-of-range values
plotData2 <- subset(plotData, grepl('^d', plotData$Detector))
#plotData2 <- droplevels(plotData2.tmp[!plotData2.tmp$Log3Concentration > 3,])
plotData2 <- droplevels(plotData2[!(plotData2$Sample=='neg' | is.na(plotData2$Sample)),])
plotData2$Detector <- factor(plotData2$Detector, levels = c('dLOC_Os03g60180', 'dLOC_Os11g06390', 'dLHS1'))

# make a data.frame for annotations
annText <- unique(plotData2[,c('Species','Detector')])

# get the positions for the labels from the data
for (i in 1:nrow(annText)){
  # get species and detector
  s <- as.character(annText$Species[i])
  d <- as.character(annText$Detector[i])
  # subset plotData2 on species and detector
  subData <- plotData2[plotData2$Species == s & plotData2$Detector == d,]
  # fit a linear model to the subset and use the gradient to make a label
  label <- coef(lm(Cq ~ Log3Concentration, data = subData))['Log3Concentration']
  annText$label[i] <- paste("gradient = ", round(label, 2), sep = '')
  # get the raw Ct values for this detector
  subData.2 <- data[substring(d,2),]
  # make labels for the minimum and maximum validated Ct values (Ct values
  # outside this range are not valid in susbsequent analysis)
  annText$labMin[i] <- round(subData.2[paste(s, '.1', sep = '')], 1)
  annText$labMax[i] <- round(subData.2[paste(s, '.', as.character(3^max(subData$Log3Concentration)), sep = '')], 1)
  # calculate positions for the labels
  annText$lPos[i] <- subset(plotData2, Detector == d & Species == s  & Log3Concentration == 0)$Cq
  annText$rPos[i] <- subset(plotData2, Detector == d & Species == s  & Log3Concentration == max(subData$Log3Concentration))$Cq
  annText$midPos[i] <- mean(annText$lPos[i], annText$rPos[i])
}

# shift the labels up where necessary
#annText[annText$Detector == 'dLOC_Os08g31980',c('lPos', 'rPos', 'midPos')] <- annText[annText$Detector == 'dLOC_Os08g31980',c('lPos', 'rPos', 'midPos')]
annText[annText$Detector == 'dLOC_Os11g06390',c('lPos', 'rPos', 'midPos')] <- annText[annText$Detector == 'dLOC_Os11g06390',c('lPos', 'rPos', 'midPos')] + 3

# re-label facets
plotData2$Detector <- factor(substr(plotData2$Detector, 2, stop = length(plotData2$Detector)), levels = c('LOC_Os03g60180', 'LOC_Os11g06390', 'LHS1'))
annText$Detector <- factor(substr(annText$Detector, 2, stop = 100), levels = c('LOC_Os03g60180', 'LOC_Os11g06390', 'LHS1'))


cairo_pdf(filename = 'ir64rpt.pdf', width = 5.059, height = 4.462, family = 'inconsolata')
g <- ggplot(data = plotData2, aes(y = Cq, x = Log3Concentration))
g + theme_minimal(base_size = 10) +
  theme(panel.border = element_rect(colour = 'grey', fill = NA)) +
  geom_point() +
  facet_grid(Detector ~ Species, scales = "free_x") +
  stat_smooth(method = lm, level = FALSE) +
  ylab(expression(Delta * C[t] * " (" * italic(vs) * ". ACT1)")) +
  xlab(expression(Log[3] * "(Concentration)")) +
  geom_text(data = annText, aes(x = 2, y = midPos, label = label), vjust = 2, size = 4) +
  geom_text(data = annText, aes(x = 0, y = lPos, label = labMin), vjust = 2, hjust = 0, size = 4) +
  geom_text(data = annText, aes(x = 4, y = rPos, label = labMax), vjust = 2, hjust = 1, size = 4)
dev.off()
ggsave('ir64rpt.pdf', height = 170*2/3, width = 257*2/4, units = 'mm', )
