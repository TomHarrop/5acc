#!/usr/bin/Rscript

library(DESeq2)
library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}
GenerateMessage("QC plots for DESeq2 results")

# set up parallel processing
cpus <- as.numeric(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
if (is.na(cpus)) {
  cpus <- 1
}
GenerateMessage(paste("Allocating", cpus, "cpu(s)"))
BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# load data
GenerateMessage("Loading base DESeq2 object")
dds.file <- "output/deseq2/dds.Rds"
if (!file.exists(dds.file)) {
  stop("Couldn't find dds.Rds")
}
dds <- readRDS(dds.file)

# pre-filter genes based on tpm cutoffs
GenerateMessage("Removing undetected genes")
detected.genes.file <- "output/tpm/detected_genes.Rds"
if (!file.exists(detected.genes.file)) {
  stop("Couldn't find detected_genes.Rds")
}
detected.genes <- readRDS(detected.genes.file)

# have to re-run deseq b/c of change in genes
dds.qc <- DESeq2::DESeqDataSetFromMatrix(
  countData = DESeq2::counts(dds)[detected.genes,],
  colData = SummarizedExperiment::colData(dds),
  design = DESeq2::design(dds))
dds.qc <- DESeq2::DESeq(dds.qc, parallel = TRUE)

vst <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
vst.assay <- SummarizedExperiment::assay(vst)

# run pca on expressed genes
pca <- prcomp(t(vst.assay))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

# set up plot
pcaPlotData <- data.frame(
  label = toupper(rownames(pca$x)),
  PCA1 = pca$x[,1],
  PCA2 = pca$x[,2],
  Stage = SummarizedExperiment::colData(vst)$stage,
  Accession = SummarizedExperiment::colData(vst)$accession
  )

library(ggplot2)
library(scales)

pcaPlot <- ggplot(pcaPlotData,
                  aes(x = PCA1, y = PCA2, colour = Accession,
                      shape = Stage, label = label)) +
  theme_grey(base_size = 8) +
  scale_color_brewer(palette = "Set1") +
  geom_point(size = 3, alpha = 0.75) +
  geom_text(colour = "black", nudge_x = 2, nudge_y = 0)
pcaPlot

# density plot
#densityPlotData <- reshape2::melt(log2(counts(dds)[exprGenes,]) + 0.5)
densityPlotData <- reshape2::melt(
  DESeq2::counts(dds.qc))
#densityPlotData <- reshape2::melt(exprVst)
densityPlotData$colour <- substr(densityPlotData$Var2, 1, 1)
densityPlot <- ggplot(densityPlotData,
       aes(x = value, fill = colour),
       alpha = 0.5) +
  xlab("Counts") + ylab("Density") + guides(fill = FALSE, colour = FALSE) +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  #scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  geom_density(alpha = 0.8, colour = NA) +
  facet_wrap(~Var2)
densityPlot

# heatmap of sample distances
sample.dists <- dist(t(vst.assay), method = "minkowski")
sample.dist.matrix <- as.matrix(sample.dists)
rownames(sample.dist.matrix) <- paste(
  rownames(SummarizedExperiment::colData(vst)),
  vst$accession, vst$stage, sep = ".")
colnames(sample.dist.matrix) <- rownames(sample.dist.matrix)


colours <- colorRampPalette(rev(RColorBrewer::brewer.pal(6, "Blues")))(255)
pheatmap::pheatmap(sample.dist.matrix,
                   clustering_distance_rows = sample.dists,
                   clustering_distance_cols = sample.dists,
                   col = colours)



