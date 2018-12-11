#!/usr/bin/env Rscript

# Table S9: DE-genes-interaction

library(data.table)

as_de_file <- "output/050_deseq/wald_tests/tfs/all/stage_accession_glaberrima.csv"
af_de_file <- "output/050_deseq/wald_tests/tfs/all/stage_accession_indica.csv"

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

fwrite(all_results, "test/Table_S9.csv")
