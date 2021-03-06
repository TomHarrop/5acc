#!/usr/bin/env Rscript

library(DESeq2)

###########
# GLOBALS #
###########


dds_file <- snakemake@input[["dds"]]
detected_genes_file <- snakemake@input[["detected_genes"]]

log_file <- snakemake@log[["log"]]
cpus <- snakemake@threads[[1]]

filtered_dds_file <- snakemake@output[["dds"]]


########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# number of CPUs to use
BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# read the dds and detected genes
detected_genes <- readRDS(detected_genes_file)
dds <- readRDS(dds_file)

# REMOVE JAPONICA FROM ALL SUBSEQUENT EXPTS
dds <- dds[, colData(dds)$accession != "japonica"]
dds$accession <- droplevels(dds$accession)

# filter
dds_filtered <- dds[detected_genes,]
design(dds_filtered) <- ~ accession + stage + accession:stage
dds_filtered <- DESeq(dds_filtered,
                      test = "LRT",
                      reduced = ~ accession + stage,
                      parallel = TRUE)

# write output
saveRDS(dds_filtered, filtered_dds_file)

# write log
sessionInfo()
