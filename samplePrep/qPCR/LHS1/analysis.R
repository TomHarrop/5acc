setwd('/home/tom/Desktop/5accessions/samplePrep/qPCR/LHS1/')

library(ReadqPCR)
library(NormqPCR)
library(ggplot2)
library(plyr)
library(grid)
library(reshape2)

source('/home/tom/functions/facetAdjust.R')

############
### MUNG ###
############

#1. read files, add species as PlateID

# easy for single plate
nipponbare <- read.table('TH 18.12.14 Nipponbare RQs.txt', skip = 1, header = TRUE, fill = TRUE, stringsAsFactors = FALSE, comment.char = '>', sep = '\t')[,1:7]
nipponbare$PlateID <- 'nipponbare'

# tricky for 2x4 plate
glaberrimaBarthii <- read.table('TH 30.1.15 Og Ob ACT1 2HK LHS1.txt', skip = 1, header = TRUE, fill = TRUE, stringsAsFactors = FALSE, comment.char = '>', sep = '\t')[,1:7]
glaberrimaBarthii$PlateID <- ''
indicaRufipogon <- read.table('TH 2.2.15 OsIR64 Or ACT1 2HK LHS1.txt', skip = 1, header = TRUE, fill = TRUE, stringsAsFactors = FALSE, comment.char = '>', sep = '\t')[,1:7]
indicaRufipogon$PlateID <- ''

# positions for each species on 2x4 plate
species1 <- c(sapply(letters[1:16], function(x) paste(toupper(x), seq(1,12), sep = '')))
species2 <- c(sapply(letters[1:16], function(x) paste(toupper(x), seq(13,24), sep = '')))

# add species names to 2x4 plate
glaberrimaBarthii[glaberrimaBarthii$Pos %in% species1,'PlateID'] <- 'glaberrima'
glaberrimaBarthii[glaberrimaBarthii$Pos %in% species2,'PlateID'] <- 'barthii'
indicaRufipogon[indicaRufipogon$Pos %in% species1,'PlateID'] <- 'indica'
indicaRufipogon[indicaRufipogon$Pos %in% species2,'PlateID'] <- 'rufipogon'

#2. Add detector names

# for single plate
LOC_Os11g06390.positions <- c(paste('A', seq(1,18), sep = ''), paste('E', seq(1,18), sep = ''))
LOC_Os03g60180.positions <- c(paste('B', seq(1,18), sep = ''), paste('F', seq(1,18), sep = ''))
ACT1.positions <- c(paste('C', seq(1,18), sep = ''), paste('G', seq(1,18), sep = ''))
LHS1.positions <- c(paste('D', seq(1,18), sep = ''), paste('H', seq(1,18), sep = ''))
nipponbare$Detector <- ''
nipponbare[nipponbare$Pos %in% LOC_Os11g06390.positions,'Detector'] <- 'LOC_Os11g06390'
nipponbare[nipponbare$Pos %in% LOC_Os03g60180.positions,'Detector'] <- 'LOC_Os03g60180'
nipponbare[nipponbare$Pos %in% ACT1.positions,'Detector'] <- 'ACT1'
nipponbare[nipponbare$Pos %in% LHS1.positions,'Detector'] <- 'LHS1'

# for 2x4 plates
act1 <- c(sapply(letters[1:16], function(x) paste(toupper(x), c(1:3, 13:15), sep = '')))
HK03 <- c(sapply(letters[1:16], function(x) paste(toupper(x), c(4:6, 16:18), sep = '')))
HK11 <- c(sapply(letters[1:16], function(x) paste(toupper(x), c(7:9, 19:21), sep = '')))
lhs1 <- c(sapply(letters[1:16], function(x) paste(toupper(x), c(10:12, 22:24), sep = '')))

glaberrimaBarthii[glaberrimaBarthii$Pos %in% act1,'Detector'] <- 'ACT1'
glaberrimaBarthii[glaberrimaBarthii$Pos %in% HK03,'Detector'] <- 'LOC_Os03g60180'
glaberrimaBarthii[glaberrimaBarthii$Pos %in% HK11,'Detector'] <- 'LOC_Os11g06390'
glaberrimaBarthii[glaberrimaBarthii$Pos %in% lhs1,'Detector'] <- 'LHS1'

indicaRufipogon[indicaRufipogon$Pos %in% act1,'Detector'] <- 'ACT1'
indicaRufipogon[indicaRufipogon$Pos %in% HK03,'Detector'] <- 'LOC_Os03g60180'
indicaRufipogon[indicaRufipogon$Pos %in% HK11,'Detector'] <- 'LOC_Os11g06390'
indicaRufipogon[indicaRufipogon$Pos %in% lhs1,'Detector'] <- 'LHS1'

#3. remove prefix from nipponbare plate

nipponbare$Name <- gsub('NIP_', '', nipponbare$Name)

#4. make a combined dataframe, tidy up

df <- rbind(nipponbare, indicaRufipogon, glaberrimaBarthii)

#remove uninteresting samples
select <- paste0(rep('N', 12), c(1:4), 'R', rep(c(1:3), each = 4))
df2 <- df[df$Name %in% select,]

#append species to sample name
df2$Name <- paste(df2$Name, df2$PlateID, sep = '.')

#5. write the munged data for reading in with readqpcr

mungedData <- data.frame(Well = df2$Pos, PlateID = df2$PlateID, Sample = df2$Name,
                         Detector = df2$Detector, Cq = df2$Cp, stringsAsFactors = FALSE)
write.table(mungedData, sep = '\t', file = 'mungedData.tab', quote = FALSE, row.names = FALSE)

###############
### ANALYSE ###
###############

data <- read.qPCR('mungedData.tab')
cts <- exprs(data)

## manually set values above maximum validated for each cutoff to Inf.

det <- rownames(cts)
samp <- colnames(cts)

#cts[grepl('ACT1', det), grepl('nipponbare', samp)][cts[grepl('ACT1', det), grepl('nipponbare', samp)] > 35] <- Inf
#cts[grepl('ACT1', det), grepl('indica', samp)][cts[grepl('ACT1', det), grepl('indica', samp)] > 39.6] <- Inf
#cts[grepl('ACT1', det), grepl('rufipogon', samp)][cts[grepl('ACT1', det), grepl('rufipogon', samp)] > 35] <- Inf
#cts[grepl('ACT1', det), grepl('glaberrima', samp)][cts[grepl('ACT1', det), grepl('glaberrima', samp)] > 34.2] <- Inf
#cts[grepl('ACT1', det), grepl('barthii', samp)][cts[grepl('ACT1', det), grepl('barthii', samp)] > 33.1] <- Inf

cts[grepl('Os03g60180', det), grepl('nipponbare', samp)][cts[grepl('Os03g60180', det), grepl('nipponbare', samp)] > 32.5] <- Inf
cts[grepl('Os03g60180', det), grepl('indica', samp)][cts[grepl('Os03g60180', det), grepl('indica', samp)] > 37.6] <- Inf
cts[grepl('Os03g60180', det), grepl('rufipogon', samp)][cts[grepl('Os03g60180', det), grepl('rufipogon', samp)] > 31.7] <- Inf
cts[grepl('Os03g60180', det), grepl('glaberrima', samp)][cts[grepl('Os03g60180', det), grepl('glaberrima', samp)] > 32] <- Inf
cts[grepl('Os03g60180', det), grepl('barthii', samp)][cts[grepl('Os03g60180', det), grepl('barthii', samp)] > 31.5] <- Inf

cts[grepl('Os11g06390', det), grepl('nipponbare', samp)][cts[grepl('Os11g06390', det), grepl('nipponbare', samp)] > 27.7] <- Inf
cts[grepl('Os11g06390', det), grepl('indica', samp)][cts[grepl('Os11g06390', det), grepl('indica', samp)] > 31.9] <- Inf
cts[grepl('Os11g06390', det), grepl('rufipogon', samp)][cts[grepl('Os11g06390', det), grepl('rufipogon', samp)] > 26.9] <- Inf
cts[grepl('Os11g06390', det), grepl('glaberrima', samp)][cts[grepl('Os11g06390', det), grepl('glaberrima', samp)] > 27.3] <- Inf
cts[grepl('Os11g06390', det), grepl('barthii', samp)][cts[grepl('Os11g06390', det), grepl('barthii', samp)] > 26.5] <- Inf

cts[grepl('LHS1', det), grepl('nipponbare', samp)][cts[grepl('LHS1', det), grepl('nipponbare', samp)] > 35] <- Inf
cts[grepl('LHS1', det), grepl('indica', samp)][cts[grepl('LHS1', det), grepl('indica', samp)] > 39.6] <- Inf
cts[grepl('LHS1', det), grepl('rufipogon', samp)][cts[grepl('LHS1', det), grepl('rufipogon', samp)] > 35] <- Inf
cts[grepl('LHS1', det), grepl('glaberrima', samp)][cts[grepl('LHS1', det), grepl('glaberrima', samp)] > 34.2] <- Inf
cts[grepl('LHS1', det), grepl('barthii', samp)][cts[grepl('LHS1', det), grepl('barthii', samp)] > 33.1] <- Inf

exprs(data) <- cts

# make case-control matrix
sampleTypes <- unique(paste(substr(colnames(cts), start = 1, stop = 2), substr(colnames(cts), start = 5, stop = 100), sep = ''))
df <- data.frame(row.names = sampleNames(data))

# this works but it's pretty nasty. searches substrings of sampleTypes in
# colnames(df) to decide whether to add a '1' or a '0'.
for (x in sampleTypes) {
  sampleProperties <- unlist(strsplit(x, split = "\\."))
  df[grepl(sampleProperties[1], rownames(df)) &
       grepl(sampleProperties[2], rownames(df)),x] <- 1
}
df[is.na(df)] <- 0
stageMatrix <- as.matrix(df)
#write.csv(stageMatrix, file = 'stageMatrix.csv') ## to visually check the matrix

# combine technical replicates
data <- combineTechReps(data)

# delta-delta-Ct with three housekeepers
hkgs <- c('LOC_Os11g06390', 'LOC_Os03g60180', 'ACT1')
#hkgs <- 'ACT1'

nipponbareSamples <- sampleTypes[grepl('nipponbare', sampleTypes)]
indicaSamples <- sampleTypes[grepl('indica', sampleTypes)]
rufipogonSamples <- sampleTypes[grepl('rufipogon', sampleTypes)]
glaberrimaSamples <- sampleTypes[grepl('glaberrima', sampleTypes)]
barthiiSamples <- sampleTypes[grepl('barthii', sampleTypes)]

nipponbareResults <- lapply(nipponbareSamples, function(x)
  deltaDeltaCq(qPCRBatch = data, hkgs = hkgs, contrastM = stageMatrix, case = x, control = "N4.nipponbare", hkgCalc = 'geom')
)
names(nipponbareResults) <- substr(nipponbareSamples, 1, 2)

indicaResults <- lapply(indicaSamples, function(x)
  deltaDeltaCq(qPCRBatch = data, hkgs = hkgs, contrastM = stageMatrix, case = x, control = "N4.indica", hkgCalc = 'geom')
)
names(indicaResults) <- substr(indicaSamples, 1, 2)

rufipogonResults <- lapply(rufipogonSamples, function(x)
  deltaDeltaCq(qPCRBatch = data, hkgs = hkgs, contrastM = stageMatrix, case = x, control = "N4.rufipogon", hkgCalc = 'geom')
)
names(rufipogonResults) <- substr(rufipogonSamples, 1, 2)

glaberrimaResults <- lapply(glaberrimaSamples, function(x)
  deltaDeltaCq(qPCRBatch = data, hkgs = hkgs, contrastM = stageMatrix, case = x, control = "N4.glaberrima", hkgCalc = 'geom')
)
names(glaberrimaResults) <- substr(glaberrimaSamples, 1, 2)

barthiiResults <- lapply(barthiiSamples, function(x)
  deltaDeltaCq(qPCRBatch = data, hkgs = hkgs, contrastM = stageMatrix, case = x, control = "N4.barthii", hkgCalc = 'geom')
)
names(barthiiResults) <- substr(barthiiSamples, 1, 2)

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

nipponbareResultsShort <- lapply(nipponbareResults, getCIs)
nipponbareFrame <- data.frame(do.call(rbind, nipponbareResultsShort), stage = names(nipponbareResultsShort), species = 'nipponbare')

indicaResultsShort <- lapply(indicaResults, getCIs)
indicaFrame <- data.frame(do.call(rbind, indicaResultsShort), stage = names(indicaResultsShort), species = 'indica')

rufipogonResultsShort <- lapply(rufipogonResults, getCIs)
rufipogonFrame <- data.frame(do.call(rbind, rufipogonResultsShort), stage = names(rufipogonResultsShort), species = 'rufipogon')

glaberrimaResultsShort <- lapply(glaberrimaResults, getCIs)
glaberrimaFrame <- data.frame(do.call(rbind, glaberrimaResultsShort), stage = names(glaberrimaResultsShort), species = 'glaberrima')

barthiiResultsShort <- lapply(barthiiResults, getCIs)
barthiiFrame <- data.frame(do.call(rbind, barthiiResultsShort), stage = names(barthiiResultsShort), species = 'barthii')

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

getPvalue <- function(resultsFrame, i) {
  return(t.testFromMeans(
    resultsFrame[i,'X2..ddCt'], # x1
    resultsFrame['N1',"X2..ddCt"], # x2
    resultsFrame[i,"X2..ddCt"] - resultsFrame[i,"X2..ddCt.min"], # sd1
    resultsFrame['N1',"X2..ddCt.max"] - resultsFrame['N1',"X2..ddCt"], # sd2
    3,3)['p'])}

nipponbareFrame$pval <- c(1, sapply(2:4, function(i) getPvalue(nipponbareFrame, i)))
indicaFrame$pval <- c(1, sapply(2:4, function(i) getPvalue(indicaFrame, i)))
rufipogonFrame$pval <- c(1, sapply(2:4, function(i) getPvalue(rufipogonFrame, i)))
glaberrimaFrame$pval <- c(1, sapply(2:4, function(i) getPvalue(glaberrimaFrame, i)))
barthiiFrame$pval <- c(1, sapply(2:4, function(i) getPvalue(barthiiFrame, i)))

# mung it all together
plotData <- rbind(nipponbareFrame, indicaFrame, rufipogonFrame, glaberrimaFrame, barthiiFrame)
rownames(plotData) <- NULL
colnames(plotData) <- c('ddCt', 'ddCt.min', 'ddCt.max', 'CI95min', 'CI95max', 'stage', 'species', 'pval')
plotData$stage <- revalue(plotData$stage, replace = c('N1' = 'Rachis\nmeristem', 'N2' = 'Branch\nmeristem',
                                                      'N3' = 'Spikelet\nmeristem', 'N4' = 'Floret\nmeristem'))
plotData$species <- revalue(plotData$species, replace = c('nipponbare' = 'O. sativa Japonica\nNipponbare', 'indica' = 'O. sativa Indica\nIR64',
                                                          'rufipogon' = 'O. rufipogon\nW1654', 'glaberrima' = 'O. glaberrima\nTog 5681',
                                                          'barthii' = 'O. barthi\nB88'))
plotData$pval[is.nan(plotData$pval)] <- 1
plotData$labels <- 3
plotData[plotData$pval < 0.05,'labels'] <- 1
plotData[plotData$pval < 0.005,'labels'] <- 2
plotData$asterisks <- c('*', '**', '')[as.numeric(plotData$labels)]

# plot 1 (lines)

g <- ggplot(plotData, aes(x = stage, y = ddCt,
                          ymin = CI95min,
                          ymax = CI95max)) +
  xlab(NULL) + ylab(expression(atop(Normalised~relative~italic("LHS1")~expression,
                                    (2^{- Delta * Delta * C[t]}*phantom(" ")%+-%" 95% CI, " *italic(n) == 3 )))) +
  stat_smooth(aes(group = species), se = FALSE, colour = 'grey') +
  geom_errorbar(width = 0.2, colour = 'grey')
g <- g + theme_minimal() + geom_point(aes(colour = stage), size = 3) +
  facet_wrap(~species, scales = 'fixed') + 
  scale_color_brewer(palette = 'Set1',
                     guide = guide_legend(label.hjust = 0, title = NULL, order = 1))

g <- g + geom_text(aes(y = CI95max, label = asterisks), vjust = 0, colour = 'grey') +
  geom_point(aes(shape = asterisks), colour = NA) +
  scale_shape_manual(values = c('1',"2",'3'),
                     labels = c(expression(paste('*'~~0.05 > italic(p), phantom() >= 0.005)),
                                expression("**"~~italic(p) < 0.005), ''),
                     guide = guide_legend(title=expression(atop(italic(t)*phantom()*"-"*test~italic(vs.),rachis~meristem~expression*":")),
                                          label.hjust = 0.5, title.hjust = 0.5)) +
  theme(legend.position = c(5/6,1/4), legend.justification = c(0.5,0.5),
        legend.key.height = unit(1.5, 'lines'),
        strip.text = element_text(face = 'bold'),
        axis.text.x = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(lineheight = 3/4))

#f <- facetAdjust(g, 'up')
cairo_pdf(filename = 'LHS1.pdf', width = 10, height = 7.5, pointsize = 12, family = "Inconsolata")
g
dev.off()

#plot 2 (bars)
b <- ggplot(plotData, aes(x = species, y = ddCt, fill = stage,
                          ymin = CI95min,
                          ymax = CI95max)) +
  theme_minimal() + theme(axis.ticks = element_blank(),
                          legend.text = element_text(lineheight = 3/4, size = 8)) +
  xlab(NULL) + ylab(expression(atop(Normalised~relative~italic("LHS1")~expression,
                                    (2^{- Delta * Delta * C[t]}*phantom(" ")%+-%" 95% CI, " *italic(n) == 3))))
b <- b + geom_bar(stat = 'identity', position = 'dodge') +
  geom_bar(stat = 'identity', position = 'dodge', colour = 'black', show_guide = FALSE) +
  scale_fill_brewer(palette = 'Set1', guide = guide_legend(title = NULL,
                                                           label.hjust = 0))
b <- b + geom_errorbar(width = 0.2, colour = 'black', position = position_dodge(0.9), alpha = 0.6) +
  geom_text(aes(y = CI95max, label = asterisks), vjust = 0, colour = 'grey', position = position_dodge(0.9))
cairo_pdf(filename = 'LHS1.bars.pdf', width = 10, height = 7.5, pointsize = 12, family = "Inconsolata")
b + scale_y_continuous(expand = c(0,0))
dev.off()

#for talk
cairo_pdf(filename = '/home/tom/Dropbox/temp/LHS1.bars.pdf',
          width = as.numeric(convertUnit(unit(218, 'mm'), unitTo = 'inches')),
          height = as.numeric(convertUnit(unit(139.25, 'mm'), unitTo = 'inches')),
          pointsize = 16, family = "Verdana")
b + scale_y_continuous(expand = c(0,0))
dev.off()

##############
### GENORM ###
##############

datagn <- data
#genorm doesn't like Inf
exprs(datagn)[exprs(datagn) == Inf] <- NA
gn <- selectHKs(datagn, method = 'geNorm', Symbols = featureNames(data), minNrHK = 2, log = FALSE)
#check that the pairwise variation is < 0.15
gn$variation

# run geNorm per-species (because this is how the normalization is done)
species <- factor(gsub('.*\\.', '', colnames(data)), levels = c('nipponbare', 'indica', 'rufipogon', 'glaberrima', 'barthii'))
gnResults <- lapply(levels(species), function(x)
  selectHKs(datagn[,species == x], method = 'geNorm', Symbols = featureNames(data), minNrHK = 2, log = FALSE))
names(gnResults) <- levels(species)

# check the ranking per species. Obviously LHS1 should lose...
sapply(levels(species), function(x) gnResults[[x]]$ranking)

# visualise the pairwise variation between 'housekeepers' (n.b. V2/3 includes
# LHS1, so not really a HK)
gnPlot <- sapply(levels(species), function(x) gnResults[[x]]$variation)
g <- ggplot(melt(gnPlot), aes(x = Var2, y = value)) + theme_minimal() +
  geom_bar(aes(fill = Var1), position = 'dodge', stat = 'identity') +
  scale_fill_brewer(palette = 'Set1') + 
  ylim(c(0,0.2)) + geom_hline(y = 0.15, linetype = 'dashed')
cairo_pdf('geNorm.pdf', width = 10/2, height = 7.5/2, pointsize = 12, family = 'Inconsolata')
g
dev.off()
