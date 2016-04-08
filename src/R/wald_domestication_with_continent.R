#!/usr/bin/env Rscript

library(DESeq2)
library(data.table)
library(ggplot2)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

# plot interesting genes
tpm <- readRDS("output/tpm/tpm_with_calls.Rds")
PlotSigGenes <- function(gene.names, counts = tpm) {
  sig.tpm <- counts[gene.names]
  sig.tpm[, gene.name := oryzr::LocToGeneName(gene)$symbols, by = gene]
  sig.tpm[is.na(gene.name), gene.name := gene]
  ggplot(sig.tpm, aes(x = stage, y = tpm, colour = species, group = species)) +
    facet_wrap(~gene.name, scales = "free_y") +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE)
}

GenerateMessage("Calculate interaction between domestication and stage")
design <- ~ continent + domestication + stage +
  continent:stage + domestication:stage +  continent:domestication +
  continent:domestication:stage

design <- ~ stage * domestication * continent

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

# pre-filter based on tpm cutoffs
GenerateMessage("Removing undetected genes")
detected.genes.file <- "output/tpm/detected_genes.Rds"
if (!file.exists(detected.genes.file)) {
  stop("Couldn't fine detected_genes.Rds")
}
detected.genes <- readRDS(detected.genes.file)
dds.domestication.with.continent <- DESeq2::DESeqDataSetFromMatrix(
  countData = DESeq2::counts(dds)[detected.genes, ],
  colData = SummarizedExperiment::colData(dds),
  design = design)

# update design
GenerateMessage("Running DESeq2 with updated design")
dds.domestication.with.continent <- DESeq2::DESeq(
  dds.domestication.with.continent, parallel = TRUE)

# extract results
resultsNames(dds.domestication.with.continent)

# stage:domestication:continent interaction
res <- results(dds.domestication.with.continent,
               name = "stageSM.domesticationdomesticated.continentAsia",
               alpha = 0.05, parallel = TRUE)

PlotSigGenes(rownames(subset(res, padj < 0.10)), counts = tpm)

oryzr::LocToGeneName(rownames(subset(res, padj < 0.01)))
res["LOC_Os03g60430",]

res.tpm <- tpm[rownames(subset(res, padj < 0.05))]
res.tpm[, gene.name := oryzr::LocToGeneName(gene)$symbols, by = gene]
res.tpm[is.na(gene.name), gene.name := gene]
library(ggplot2)
ggplot(res.tpm, aes(x = stage, y = tpm, colour = species, group = species)) +
  facet_wrap(~gene.name, scales = "free_y") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

# domestication effect
dom <- results(dds.domestication.with.continent,
               contrast = list(
                 c("domestication_domesticated_vs_wild",
                   "stageSM.domesticationdomesticated",
                   "domesticationdomesticated.continentAsia")
               ), alpha = 0.05, parallel = TRUE)
oryzr::LocToGeneName(rownames(head(dom[order(dom$padj),], 20)))
PlotSigGenes(rownames(head(dom[order(dom$padj),], 30)), counts = tpm)

# stage effect
st <- results(dds.domestication.with.continent,
               contrast = list(
                 c("stage_SM_vs_PBM",
                   "stageSM.domesticationdomesticated",
                   "stageSM.continentAsia")
               ), alpha = 0.05, parallel = TRUE)
PlotSigGenes(rownames(head(st[order(st$padj),], 30)), counts = tpm)
summary(st)

# stage:domestication in Africa
sd.africa <- results(dds.domestication.with.continent,
                     name = "stageSM.domesticationdomesticated",
                     alpha = 0.05, parallel = TRUE)

# stage:domestication in Asia (interaction term plus second-order interaction)
sd.asia <- results(dds.domestication.with.continent,
                   contrast = list(
                     c("stageSM.domesticationdomesticated",
                       "stageSM.domesticationdomesticated.continentAsia")),
                   alpha = 0.05, parallel = TRUE)

# genes in both african and asian lists
sd.africa.table <- data.table(data.frame(sd.africa), keep.rownames = TRUE, key = "rn")
sd.asia.table <- data.table(data.frame(sd.asia), keep.rownames = TRUE, key = "rn")
sd.combined <- sd.africa.table[sd.asia.table, .(
  gene.id = rn,
  l2fc.africa = log2FoldChange,
  l2fc.asia = i.log2FoldChange,
  padj.africa = padj,
  padj.asia = i.padj
)]
sd.combined[padj.africa < 0.05]
sd.combined[padj.asia < 0.05]
sd.combined[padj.africa < 0.1 & padj.asia < 0.1, PlotSigGenes(gene.id)]

# what is the proper interaction term for stage:domestication???
attr(dds.domestication.with.continent, "modelMatrix")
sum.all <- results(dds.domestication.with.continent,
              contrast = list(
                c("stageSM.domesticationdomesticated"),
                c("stageSM.domesticationdomesticated.continentAsia")),
              alpha = 0.05, parallel = TRUE)
oryzr::LocToGeneName(rownames(subset(sum.all, padj < 0.05)))
oryzr::LocToGeneName(rownames(subset(r2, padj < 0.05)))$symbols

tpm <- readRDS('output/tpm/tpm_with_calls.Rds')
r2.tpm <- tpm[rownames(subset(r2, padj < 0.05))]
r2.tpm[, gene.name := oryzr::LocToGeneName(gene)$symbols, by = gene]
r2.tpm[is.na(gene.name), gene.name := gene]
library(ggplot2)
ggplot(r2.tpm, aes(x = stage, y = tpm, colour = species, group = species)) +
  facet_wrap(~gene.name, scales = "free_y") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)



summary(r2)
plotCounts(dds.domestication.with.continent, "LOC_Os01g17396", intgroup = c("accession", "stage"))

summary(res)
bof <- head(res[order(res$padj), ])
oryzr::LocToGeneName(rownames(bof))
data.table(oryzr::LocToGeneName(rownames(subset(res, padj < 0.05))),keep.rownames = TRUE)

plotCounts(dds.domestication.with.continent, "LOC_Os03g60430",
           intgroup = c("accession", "stage"))

