#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

# Table S8: clustered-genes

cluster_file <- snakemake@input[["clusters"]]

# read data
clusters <- fread(cluster_file)

clusters[, c("RapID",
             "names",
             "OgroObjective",
             "OgroRef") := NULL]
setnames(clusters, "MsuID", "gene_id")
setcolorder(clusters, c("cluster", "gene_id"))

fwrite(clusters, snakemake@output[["table1"]])

sessionInfo()