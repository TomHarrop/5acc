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
getSpecRes <- function(spec) {
  message(paste("Generating results for", spec))
  contrast <- list(
    c(paste0("stageN3.species", spec)),
    c(paste0("stageN2.species", spec)))
  return(DESeq2::results(dds, contrast))
}
speciesResults <- lapply(sampleToSpecies, getSpecRes)
names(speciesResults) <- sampleToSpecies

# quick explore
plotData <- data.frame(
  row.names = rownames(speciesResults[[1]]),
  sapply(speciesResults, function(x) x$log2FoldChange))

plotData <- plotData[complete.cases(plotData),]

comparisons <- list(
  c('barthii', 'glaberrima'),
  c('rufipogon', 'Japonica'),
  c('rufipogon', 'Indica'))

getCompLfcs <- function(comparison) {
  return(data.frame(
    gene_name = rownames(plotData),
    comparison = paste0(comparison[1], 'vs', comparison[2]),
    s1.name = comparison[1],
    s2.name = comparison[2],
    s1.lfc = plotData[,comparison[1]],
    s2.lfc = plotData[,comparison[2]]))
}

compData <- lapply(comparisons, getCompLfcs)
compData <- do.call(rbind, compData)

ggplot(compData, aes(x = s1.lfc, y = s2.lfc)) +
  geom_point() +
  facet_wrap(~comparison)

# filter on variance
varOrder <- apply(plotData, 1, var)
plotData.filter <- plotData[rev(order(varOrder)),][c(1:1000),]

# set up the expressionSet
pData <- data.frame(
  species = as.factor(colnames(plotData)),
  row.names = colnames(plotData))

phenoData <- new('AnnotatedDataFrame', data = pData)
mf.e <- ExpressionSet(assayData = as.matrix(plotData.filter), phenoData = phenoData)

# run clustering
mf.std <- Mfuzz::standardise(mf.e)
m1 <- Mfuzz::mestimate(mf.std)
Mfuzz::Dmin(mf.std, m = m1, crange = seq(2, 20, 1), repeats = 3, visu = TRUE)
set.seed(1)
c1 <- mfuzz(mf.std, centers = 9, m = m1)

# get max membership per gene
alphaCores <- acore(mf.std, c1, min.acore = 0.5)
maxMems.list <- sapply(rownames(mf.std), function(gid)
  as.numeric(which.max(sapply(alphaCores, function(x) x[gid,'MEM.SHIP']))))
cluster <- as.data.table(unlist(maxMems.list), keep.rownames = TRUE)

# make plots
extrafont::loadfonts()

pd.dt <- as.data.table(compData)
setkey(pd.dt, 'gene_name')
setkey(cluster, V1)

pd.dt.cluster <- cluster[pd.dt]
pd.dt.cluster[, gene_id := V1][, cluster := V2][,c('V1', 'V2') := NULL]

pd.dt.cluster[,comparison := plyr::revalue(comparison,
              replace = c(
                "barthiivsglaberrima" = expression(italic(O.~barthii~vs.~O.~glaberrima)),
                "rufipogonvsJaponica" = expression(italic(O.~rufipogon~vs.~O.~sativa)*" ssp. "*italic(japonica)),
               "rufipogonvsIndica" = expression(italic(O.~rufipogon~vs.~O.~sativa)*" ssp. "*italic(indica))
              ))]

g <- ggplot() +
  geom_point(data = pd.dt.cluster[is.na(cluster)], aes(x = s1.lfc, y = s2.lfc),
             colour = 'grey', alpha = 0.2, size = 5) +
  geom_point(data = pd.dt.cluster[!is.na(cluster)],
             aes(x = s1.lfc, y = s2.lfc, colour = as.factor(cluster)),
             alpha = 0.8, size = 5) +
  facet_grid(.~comparison, labeller = label_parsed) +
  scale_colour_brewer(palette = 'Set1', guide = FALSE) + 
  xlab(expression(log[2](Fold~Change)*" in wild species")) +
  ylab(expression(log[2](Fold~Change)*" in domesticated species")) +
  theme_minimal(base_size = 24, base_family = 'Lato') +
  coord_fixed() +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 10),
        rect = element_rect(fill = 'transparent', colour = NA))

# output
pdfLocation <- paste0('fig/', outputBasename, '.pdf')
logLocation <- paste0('log/', outputBasename, '.sessionInfo.txt')

wA <- grid::convertUnit(grid::unit(800, 'mm'), unitTo = 'in', valueOnly = TRUE)
hA <- grid::convertUnit(grid::unit(400, 'mm'), unitTo = 'in', valueOnly = TRUE)
cairo_pdf(pdfLocation, width = wA, height = hA,
          family = 'Lato', bg = 'transparent')
g
dev.off()

# logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
  paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
  capture.output(sessionInfo()))

writeLines(sInf, logLocation)
