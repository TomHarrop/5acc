#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# Table S6: DE-genes-stages

library(data.table)

de_genes_file <- snakemake@input[["de_genes"]]

de_genes <- fread(de_genes_file)
de_genes[, c("design",
             "baseMean",
             "stat",
             "wald_test",
             "RapID",
             "names",
             "OgroObjective",
             "OgroRef") := NULL]
setorder(de_genes, padj, na.last = TRUE)

fwrite(de_genes, snakemake@output[["table1"]])

sessionInfo()
