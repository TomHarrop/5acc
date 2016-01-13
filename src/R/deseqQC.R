#!/usr/bin/Rscript

library(DESeq2)

# load the deseq files
ddsFile <- "output/deseq2/dds.Rds"
vstFile <- "output/deseq2/vst.Rds"

findFiles <- function(x) {
  if (!file.exists(x)) {
    stop("Couldn't find dds.Rds, exiting\n")
    quit(save = "no", status = 1)
  }
}
lapply(list(ddsFile, vstFile), findFiles)

dds <- readRDS(ddsFile)
vst <- readRDS(vstFile)

# choose genes for vst; come back to this once expression cutoffs have been defined
qS <- quantile(rowSums(counts(dds)), 0.7)
qM <- quantile(rowMeans(counts(dds)), 0.7)

exprGenes <- rownames(counts(dds)[rowSums(counts(dds)) > qS | rowMeans(counts(dds)) > qM, ])
exprVst <- assay(vst[exprGenes,])

# run pca on expressed genes
pca <- prcomp(t(exprVst))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

# set up plot
pcaPlotData <- data.frame(
  label = toupper(rownames(pca$x)),
  PCA1 = pca$x[,1],
  PCA2 = pca$x[,2],
  Stage = GenomicRanges::colData(vst)$stage,
  Accession = GenomicRanges::colData(vst)$accession
  )

library(ggplot2)
library(scales)

pcaPlot <- ggplot(pcaPlotData,
                  aes(x = PCA1, y = PCA2, colour = Accession,
                      shape = Stage, label = label)) +
  theme_grey(base_size = 8) +
  scale_color_brewer(palette = "Set1") +
  #geom_point(size = 5, alpha = 0.8) +
  #geom_label(size = 2.82, colour = "black", nudge_x = 2, nudge_y = 0)
  geom_text(size = 2.82)
pcaPlot

# density plot
#densityPlotData <- reshape2::melt(log2(counts(dds)[exprGenes,]) + 0.5)
densityPlotData <- reshape2::melt(counts(dds)[exprGenes,])
#densityPlotData <- reshape2::melt(exprVst)
densityPlotData$colour <- substr(densityPlotData$Var2, 1, 1)
densityPlot <- ggplot(densityPlotData,
       aes(x = value, fill = colour),
       alpha = 0.5) +
  xlab("Transformed counts") + ylab(NULL) + guides(fill = FALSE, colour = FALSE) +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  #scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  geom_density(alpha = 0.8, colour = NA) +
  facet_wrap(~Var2)

