#!/usr/bin/env Rscript

library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

GenerateMessage("Generate tables of wald test results")

# data files
list.files("output/deseq2/wald_tests",
           pattern = ".*results_table.Rds")

dom.asia <- readRDS(
  "output/deseq2/wald_tests/domestication_asia_results_table.Rds")
dom.continent <- readRDS(
  "output/deseq2/wald_tests/domestication_by_continent_results_table.Rds")

head(dom.asia)
head(dom.continent)

dom.all <- rbind(dom.asia, dom.continent)

dom.all[padj < 0.05 & domestication %in% c("japonica", "indica")]


# wide (p by domestication)
p.by.dom <- dcast(dom.all, gene ~ domestication, value.var = "padj")
p.by.dom[indica < 0.05 & japonica < 0.05, oryzr::LocToGeneName(gene)]
p.by.dom[africa < 0.05 & (indica < 0.05 | japonica < 0.05), oryzr::LocToGeneName(gene)]
p.by.dom["LOC_Os03g60430"]

# ignore P, take genes were lfc is above 1.5
dom.all[, min.l2fc := abs(log2FoldChange) - lfcSE]
ml2fc.by.dom <- dcast(dom.all, gene ~ domestication, value.var = "min.l2fc")
ml2fc.by.dom[indica >= log(1.5, 2) &
               japonica >= log(1.5, 2) &
               africa > log(1.5, 2),
             oryzr::LocToGeneName(gene)]

goi <- ml2fc.by.dom[indica >= log(1.5, 2) &
               japonica >= log(1.5, 2) &
               africa > log(1.5, 2),
             unique(gene)]
