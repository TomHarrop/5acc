#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# Table S9: DE-genes-interaction

library(data.table)

as_de_file <- snakemake@input[["as_de"]]
af_de_file <- snakemake@input[["af_de"]]

all_results <- rbindlist(list(Asia = fread(as_de_file),
                              Africa = fread(af_de_file)),
                         idcol = "Continent")
all_results[, c("design",
                "wald_test",
                "baseMean",
                "stat",
                "RapID",
                "names",
                "OgroObjective",
                "OgroRef") := NULL]
all_results[, mpv := mean(padj), by = gene_id]
setorder(all_results, mpv, Continent, na.last = TRUE)
all_results[, mpv := NULL]

fwrite(all_results, snakemake@output[["table1"]])
