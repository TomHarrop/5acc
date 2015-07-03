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

# find the most recent cutadapt output
outputDirs <- list.dirs(path = 'output', full.names = TRUE, recursive = FALSE)
cutadaptDir <- rev(sort(outputDirs[grep('cutadapt', outputDirs)]))[1]

# find the most recent STAR output
outputDirs <- dir(path = cutadaptDir, pattern = "STAR", full.names = TRUE)
starDir <- rev(sort(outputDirs[grep('STAR', outputDirs)]))[1]

# get the read count files
fileName <- dir(starDir, pattern = "ReadsPerGene.out.tab", full.names = TRUE)

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

# generate countData
names(fileName) <- sampleName
countData <- sapply(fileName, function(x)
  read.delim(x, header = FALSE, stringsAsFactors = FALSE)[[4]])
rownames(countData) <- read.delim(fileName[[1]], header = FALSE,
                                  stringsAsFactors = FALSE)[[1]]
# generate colData
colData <- data.frame(row.names = sampleName, stage, species)

# read matrix into DESeq2, excluding the rows starting with "N_"
dds <- DESeq2::DESeqDataSetFromMatrix(countData = subset(countData, !grepl("^N_", rownames(countData))),
                                      colData = colData,
                                      design = ~ stage + species + stage:species)

# run DESeq2
dds <- DESeq2::DESeq(dds)

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
c1 <- mfuzz(mf.std, centers = 7, m = m1)

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
logLocation <- paste0('log/', outputBasename, '.sessionInfo.txt')

ggsave(pdfLocation, width = 16, height = 8,
       units = 'in')

writeLines(capture.output(sessionInfo()), logLocation)

# not run
head(
plotData[rev(order(rowSums(plotData))),]
)
DESeq2::plotCounts(dds, 'LOC_Os08g29854', intgroup = c('species', 'stage'), transform = TRUE)

