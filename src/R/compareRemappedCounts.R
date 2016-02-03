#!/usr/bin/Rscript

# 1. DEPENDENCIES

library(data.table)
library(DESeq2)
library(ggplot2)

# 2. FUNCTIONS

TidyCountData <- function(x, count.data.list) {
  # Selects column V4 (stranded = reverse) and returns tidy count data for input
  # into DESeq2.
  #
  # Args:
  #   x (chr): a name in count.data.raw
  #   count.data.list (list): list of data.frames read from quant.files
  clean.counts <- data.table::copy(count.data.list[[x]])
  rownames(clean.counts) <- clean.counts$V1
  data.table::setnames(clean.counts, old = "V4", new = x)
  return(subset(clean.counts, select = x))
}

# 3. CODE

# find all the quant files
quant.files.all <- list.files("output", pattern = "ReadsPerGene", recursive = TRUE,
                         full.names = TRUE)

# load the normal mappings
quant.files <- quant.files.all[!grepl("remap", quant.files.all)]
if (length(quant.files) == 0) {
  stop("Couldn't find quant files, exiting\n")
  quit(save = "no", status = 1)
}
count.data.raw <- lapply(quant.files, read.table, stringsAsFactors = FALSE)

# make a table of countData
names(count.data.raw) <- gsub(".*([[:upper:]][[:digit:]]).*", "\\1", quant.files)
tidy.data.list <- lapply(names(count.data.raw), TidyCountData,
                         count.data.list = count.data.raw)
count.data <- do.call(cbind, tidy.data.list)

# load the remapped counts
remapped.quant.files <- quant.files.all[grepl("remap", quant.files.all)]
remapped.counts.raw <- lapply(remapped.quant.files, read.table,
                              stringsAsFactors = FALSE)
names(remapped.counts.raw) <- gsub(".*([[:upper:]][[:digit:]]).*", "\\1",
                                   remapped.quant.files)
tidy.remapped.counts <- lapply(names(remapped.counts.raw), TidyCountData,
                               count.data.list = remapped.counts.raw)
remapped.count.data <- do.call(cbind, tidy.remapped.counts)

# make a data.frame of combined counts
count.data.table <- data.table(count.data, keep.rownames = TRUE)
remapped.count.data.table <- data.table(remapped.count.data, keep.rownames = TRUE)
count.data.table.combined <- rbindlist(list(count.data.table,
                                            remapped.count.data.table), fill = TRUE)[,
                                                    lapply(.SD, sum, na.rm = TRUE),
                                                    by = rn]
count.data.combined <- data.frame(count.data.table.combined, row.names = "rn")
setnames(count.data.combined, names(count.data.combined), 
         paste(names(count.data.combined), "combined", sep = "."))

# join combined counts with original counts
count.data.full <- cbind(count.data, count.data.combined)

# make colData for DESeq2 object
col.data.full <- data.table(
  rn = names(count.data.full)
)
col.data.full[, accession := plyr::mapvalues(substr(rn, 1, 1),
                                             from = c("J", "I", "R", "G", "B"),
                                             to = c("japonica", "indica", "rufipogon",
                                                    "glaberrima", "barthii"))]
col.data.full[as.numeric(substr(rn, 2, 2)) < 4, stage := "PBM"]
col.data.full[as.numeric(substr(rn, 2, 2)) > 3, stage := "SM"]
col.data.full[grep("combined", rn), combined := "yes"]
col.data.full[is.na(combined), combined := "no"]
col.data.full.df <- data.frame(col.data.full, row.names = "rn")
ind <- sapply(col.data.full.df, is.character)
col.data.full.df[ind] <- lapply(col.data.full.df[ind], factor)



# make a DESeqDataSet from the full count table
dds.full <- DESeqDataSetFromMatrix(
  countData = count.data.full[!grepl("^N", rownames(count.data.full)),],
  colData = col.data.full.df,
  design = ~ stage + accession + combined + combined:stage
)

# run DESeq2 and extract results for "combined" factor
dds.full <- DESeq(dds.full)
resultsNames(dds.full)
res <- results(dds.full, name = "combined_yes_vs_no")

# set up normalized counts for plotting
counts.normalized.wide <- counts(dds.full, normalized = TRUE)
counts.normalized <- data.table(melt(counts.normalized.wide))
setnames(counts.normalized, names(counts.normalized),
         c("msu.id", "library", "counts"))
setkey(counts.normalized, "library")
setkey(col.data.full, "rn")
plot.data <- col.data.full[counts.normalized]

PlotGeneCounts <- function(gene.name, pd = plot.data) {
  pd <- plot.data
  gene.plot.data <- pd[msu.id == gene.name]
  g <- ggplot(data = gene.plot.data,
              mapping = aes(x = stage, y = counts, group = accession)) +
    geom_smooth(method = lm, se = FALSE, na.rm = TRUE, size = 0.5) +
    geom_point(position = position_jitter(width = 0.2), size = 4, alpha = 0.5) +
    facet_grid(accession ~ combined)
  return(g)  
}

PlotGeneCounts("LOC_Os07g41970")
oryzr::LocToGeneName("LOC_Os07g41970")

summary(res)

subset(res, padj < 0.05)
res[order(abs(res$log2FoldChange), decreasing = TRUE), ]

counts(dds.full)["LOC_Os03g63250", ]
counts(dds.full)["LOC_Os03g63250", grepl("^G", colnames(dds.full))]

oryzr::LocToGeneName("LOC_Os03g63250")

plotCounts(dds.full, "LOC_Os12g34620", intgroup = c("accession", "stage", "combined"))

