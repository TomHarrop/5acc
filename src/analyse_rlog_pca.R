# this script:
# - Selects genes that do pass DESeq2 quality cutoff
# - rlog transform them
# - takes the PCA

# and saves both rlog and PCA


# Set up environment ------------------------------------------------------

library(tidyverse)
library(DESeq2)

# Load deseqdataset and calculate rlog ------------------------------------

dds <- readRDS("../data-raw/dds.Rds")
design(dds) <- ~ accession + stage + accession:stage
dds <- DESeq(dds,
             test = "LRT",
             reduced = ~ accession + stage,
             parallel = TRUE)
good_ones <- !is.na(results(dds)$padj)
sum(!good_ones) # Deseq removes 3438 genes for quality reasons (low counts?????)

# Do I have to take the rlog again after I clean out the bad ones?
dds <- dds[good_ones, ]
rld <- rlog(dds, blind = FALSE)

# PCA on rlog counts ------------------------------------------------------

# pca on transposed dataset, as they do it in DESeq2

rlog_counts <- assay(rld)

pc <- prcomp(t(rlog_counts ), center = TRUE)
pcx <- as.data.frame(pc$x)
pcx$ID <- rownames(pcx)
pcx <- as_tibble(pcx)

# Sample description for plotting PCA
samples <- colData(rld)
samples$ID <- rownames(samples)
samples <- as_tibble(as.data.frame(samples))
pcx <- inner_join(samples, pcx, by = "ID")

# save rotations for set enrichments
pcro <- pc$rotation
pcro <- sweep(pcro, 2, colSums(pcro), "/")
pcro <- as.data.frame(pcro); pcro$locus_id <- rownames(pcro)


# save everything ---------------------------------------------------------

save(rld, pc, pcx, pcro, file = "../data/rlog-pca.Rdata")



