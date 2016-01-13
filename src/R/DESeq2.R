#!/usr/bin/Rscript

library(data.table)
library(DESeq2)
library(ggplot2)

# load the quant files
quantFiles <- list.files("output", pattern = "ReadsPerGene", recursive = TRUE,
                         full.names = TRUE)
if (length(quantFiles) == 0) {
  stop("Couldn't find quant files, exiting\n")
  quit(save = "no", status = 1)
}
rawCountData <- lapply(quantFiles, read.table, stringsAsFactors = FALSE)

# make a table of countData for DESeq2
names(rawCountData) <- gsub(".*([[:upper:]][[:digit:]]).*", "\\1", quantFiles)
tidyCountData <- function(x) {
  cleanCounts <- data.table::copy(rawCountData[[x]])
  rownames(cleanCounts) <- cleanCounts$V1
  data.table::setnames(cleanCounts, old = "V4", new = x)
  return(subset(cleanCounts, select = x))
}
countDataList <- lapply(names(rawCountData), tidyCountData)
countData <- do.call(cbind, countDataList)

# set up colData
colData.table <- data.table(rn = colnames(countData))
colData.table[, accession := plyr::mapvalues(substr(rn, 1, 1),
                              from = c("J", "I", "R", "G", "B"),
                              to = c("japonica", "indica", "rufipogon",
                                     "glaberrima", "barthii"))]
colData.table[as.numeric(substr(rn, 2, 2)) < 4, stage := "PBM"]
colData.table[as.numeric(substr(rn, 2, 2)) > 3, stage := "SM"]
colData.table[accession %in% c("barthii", "rufipogon"), domestication := "wild"]
colData.table[!accession %in% c("barthii", "rufipogon"), domestication := "domesticated"]

colData <- data.frame(colData.table, row.names = "rn")
ind <- sapply(colData, is.character)
colData[ind] <- lapply(colData[ind], factor)
colData$domestication <- relevel(colData$domestication, "wild")
colData$accession <- relevel(colData$accession, "japonica")

# build DESeq2 object
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = subset(countData, !grepl("^N", rownames(countData))),
  colData = colData,
  design = ~ domestication + stage + domestication:stage)
# separate object for comparisons between species
ddsSpecies <- DESeq2::DESeqDataSetFromMatrix(
  countData = subset(countData, !grepl("^N", rownames(countData))),
  colData = colData,
  design = ~ accession + stage + accession:stage)

# run DESeq2
dds <- DESeq2::DESeq(dds)
ddsSpecies <- DESeq2::DESeq(ddsSpecies)

# investigate stage:domestication interaction
DESeq2::resultsNames(dds)
res <- DESeq2::results(dds, name = "domesticationdomesticated.stageSM", alpha = 0.05)
summary(res)
someGenes <- rownames(subset(as.data.frame(res), padj < 0.05))

oryzr::LocToGeneName(someGenes)

DESeq2::plotCounts(dds, gene = "LOC_Os05g03760", intgroup = c("stage", "accession"))
DESeq2::plotCounts(dds, gene = "LOC_Os04g11980", intgroup = c("stage", "accession"))

# investigate species-specific changes
DESeq2::resultsNames(ddsSpecies)
someRes <- DESeq2::results(ddsSpecies, contrast = c("accession", "japonica", "rufipogon"), lfcThreshold = 5, alpha = 0.05)
subset(data.frame(someRes), padj < 0.05)
rufipogon <- DESeq2::results(ddsSpecies, name = "accessionrufipogon.stageSM", alpha = 0.05, lfcThreshold = 1)
summary(rufipogon)
oryzr::LocToGeneName(rownames(subset(rufipogon, padj < 0.05)))


oryzr::LocToGeneName(someMoreGenes)

# some plots
normCounts.wide <- data.frame(DESeq2::counts(dds, normalized = TRUE))
normCounts.wide$msuId <- rownames(normCounts.wide)
normCounts <- reshape2::melt(normCounts.wide, id.vars = "msuId",
                             variable.name = "Library", value.name = "count")
normCounts <- as.data.table(normCounts, keep.rownames = F)
setkey(normCounts, "Library")
setkey(colData.table, "rn")

plotData <- colData.table[normCounts]
setnames(plotData, "rn", "Library")

ggplot(plotData[msuId == "LOC_Os11g10590"], aes(x = stage, y = count)) +
  stat_smooth(aes(group = accession), method = "lm", se = FALSE, size = 0.5) +
  geom_point(position = position_jitter(width = 0.2)) +
  facet_wrap(~ accession)
oryzr::LocToGeneName("LOC_Os02g47390")
 