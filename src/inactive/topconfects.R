#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(DESeq2)
library(topconfects)

dds_file <- "output/050_deseq/filtered_dds.Rds"
cpus <- 8

BiocParallel::register(BiocParallel::MulticoreParam(cpus))

dds <- readRDS(dds_file)
design(dds) <- ~ stage + accession
dds <- DESeq(dds)
dconfects <- deseq2_confects(dds,
                             name  = "stage_SM_vs_PBM",
                             step = 0.01,
                             fdr = 0.01)
dconf_dt <- data.table(dconfects$table)

dconf_dt[abs(confect) < log(1.5, 2)]
dconf_dt[!is.na(confect) & abs(confect) > log(1.5, 2), unique(name)]

