#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)

###########
# GLOBALS #
###########

tfdb_file <- snakemake@input[["tfdb"]]
det_genes_file <- snakemake@input[["detected_genes"]]
dds_file <- snakemake@input[["dds"]]

log_file <- snakemake@log[["log"]]

dds_tfs_only <- snakemake@output[["dds"]]

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read the dds
dds <- readRDS(dds_file)

# fix factors, thanks a lot DESeq
dds$continent <- factor(dds$continent)
dds$domestication <- factor(dds$domestication)

# get the list of detected genes
det_genes <- readRDS(det_genes_file)

# get a list of transcription factors
tfdb <- readRDS(tfdb_file)
all_tfs <- tfdb[, unique(`Protein ID`)]

# subset the dds for expressed TFs
kept_genes <- intersect(names(dds), intersect(all_tfs, det_genes))
dds_exp_tfs <- dds[kept_genes]

# REMOVE JAPONICA FROM ALL SUBSEQUENT EXPTS
dds_exp_tfs <- dds_exp_tfs[, colData(dds_exp_tfs)$accession != "japonica"]
dds_exp_tfs$accession <- droplevels(dds_exp_tfs$accession)

# write output
saveRDS(dds_exp_tfs, dds_tfs_only)

# write log
sessionInfo()
