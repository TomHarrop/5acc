#!/usr/bin/Rscript

library(ggplot2)
library(Mfuzz)
library(data.table)

# set variables
scriptName <- 'quickFig'
outputBasename <- paste(
  Sys.Date(),
  scriptName,
  sep = "-"
)


# get a list of htseq files
htdirs <- dir("output/", pattern = 'htseq')
directory <- paste0('output/', htdirs[length(htdirs)])
fileName <- dir(directory, pattern = ".*htseq")

# generate sampleName
sampleName <- gsub('.*([A-Z][0-9]).*', '\\1', fileName)

# generate stage column
sampleToStage <- rep(c('N2', 'N3'), each = 3)
stage <- sampleToStage[as.numeric(gsub('[^0-9]', '', sampleName))]

# generate species column
sampleToSpecies <- c(
  B = "barthii",
  G = "glaberrima",
  J = "Japonica",
  I = "Indica",
  R = "rufipogon")
species <- sampleToSpecies[gsub('[^BGJIR]', '', sampleName)]

# make sample table
sampleTable <- data.frame(sampleName, fileName, stage, species)

# import htseq results
dds <- DESeq2::DESeqDataSetFromHTSeqCount(sampleTable, directory, design = ~ stage + species + stage:species, )

# run DESeq2
dds <- DESeq2::DESeq(dds)
saveRDS(dds, 'tmp/dds.Rds')
dds <- readRDS('tmp/dds.Rds')

# extract results
contrast <- list(
  c('stageN3.speciesglaberrima', 'stageN2.speciesglaberrima'),
  c('stageN3.speciesbarthii', 'stageN2.speciesbarthii')    
)
africa <- DESeq2::results(dds, contrast)

contrast <- list(
  c('stageN3.speciesJaponica', 'stageN2.speciesJaponica'),
  c('stageN3.speciesrufipogon', 'stageN2.speciesrufipogon')    
)
japonica <- DESeq2::results(dds, contrast)

contrast <- list(
  c('stageN3.speciesIndica', 'stageN2.speciesIndica'),
  c('stageN3.speciesrufipogon', 'stageN2.speciesrufipogon')  
)
indica <- DESeq2::results(dds, contrast)

# quick explore
plotData <- data.frame(
  row.names = rownames(africa),
  africa = africa$log2FoldChange,
  japonica = japonica$log2FoldChange,
  indica = indica$log2FoldChange
)

plotData <- plotData[complete.cases(plotData),]

plotData.melt <- reshape2::melt(plotData, variable.name = 'comparison', value.name = 'l2fc')

ggplot(plotData.melt, aes(x = l2fc, group = comparison, fill = comparison)) +
  xlim(-0.0001, 0.0001) +
  geom_density(alpha = 0.2)

# filter on variance
varOrder <- apply(plotData, 1, var)
plotData.filter <- plotData[rev(order(varOrder)),][c(1:1000),]

# set up the expressionSet
pData <- data.frame(
  comparison = as.factor(c('africa', 'japonica', 'indica')),
  row.names = c('africa', 'japonica', 'indica')
)
phenoData <- new('AnnotatedDataFrame', data = pData)
mf.e <- ExpressionSet(assayData = as.matrix(plotData.filter), phenoData = phenoData)

# run some clustering
mf.std <- standardise(mf.e)
m1 <- mestimate(mf.std)
Dmin(mf.std, m = m1, crange = seq(2, 9, 1), repeats = 1, visu = TRUE)
set.seed(1)
c1 <- mfuzz(mf.std, centers = 6, m = m1)

# make some plots

pd.dt <- as.data.table(plotData, keep.rownames = TRUE)
pd.dt[, gene_id := rn][, rn := NULL]
setkey(pd.dt, 'gene_id')

cluster <- as.data.table(c1$cluster, keep.rownames = TRUE)
setkey(cluster, V1)

pd.dt.cluster <- cluster[pd.dt]
pd.dt.cluster[, gene_id := V1][, cluster := V2][,c('V1', 'V2') := NULL]

pd.melt <- reshape2::melt(pd.dt.cluster,
                          id.vars = c('gene_id', 'africa', 'cluster'),
                          measure.vars = c('japonica', 'indica'))

ggplot() +
  geom_point(data = pd.melt[is.na(cluster)],
               aes(y = value, x = africa), colour = 'grey', alpha = 0.4) +
  geom_point(data = pd.melt[!is.na(cluster)],
             aes(y = value, x = africa, colour = as.factor(cluster)), alpha = 0.6) +
  facet_wrap(~ variable) +
  scale_colour_brewer(palette = 'Set1', guide = FALSE) +
  xlab(expression("Relative change in expression"~(italic(O.~glaberrima~vs.~O.barthii)))) +
  ylab(expression("Relative change in expression"~(italic(O.~sativa~vs.~O.rufipogon)))) +
  theme_minimal(base_size = 16) +
  theme(strip.text = element_text(face = 'italic'),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 10))

# output

pdfLocation <- paste0('fig/', outputBasename, '.pdf')
logLocation <- paste0('fig/', outputBasename, '.sessionInfo.txt')

ggsave(pdfLocation, width = 16, height = 8,
       units = 'in')

writeLines(capture.output(sessionInfo()), logLocation)
