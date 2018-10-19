#!/usr/bin/env Rscript

# Table S6: DE-genes-stages

library(data.table)

de_genes_file <- "output/050_deseq/wald_tests/expr_genes/all/stage.csv"

de_genes <- fread(de_genes_file)
de_genes[, c("design", "wald_test", "RapID", "names", "OgroObjective", "OgroRef") := NULL]
setorder(de_genes, padj, na.last = TRUE)

fwrite(de_genes, "test/Table_S6.csv")
