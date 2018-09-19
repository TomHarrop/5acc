#!/usr/bin/Rscript

library(data.table)
library(ggplot2)

source("src/R/column_plots.R")

# get tpm values
tpm.file <- "output/tpm/tpm_with_calls.Rds"
if (!file.exists(tpm.file)) {
  stop("Couldn't find tpm.Rds")
}
tpm.long <- readRDS(tpm.file)

ord <- c("O. rufipogon", "O. sativa indica", "O. sativa japonica",
         "O. barthii", "O. glaberrima")
tpm.long[, accession := substr(sample, 1, 1)]
tpm.long[ , accession := factor(plyr::mapvalues(
  accession,
  from = c("R", "I", "J", "B", "G"),
  to = ord),  levels = ord)
  ]

hirano.file <- "data/goi/hirano_hormone_genes.Rds"
hirano.genes <- readRDS(hirano.file)

hirano.genes[, unique(hormone.class)]
goi <- hirano.genes[hormone.class == "ABA" & !is.na(msu.id), unique(msu.id)]

ColumnPlot(goi)

# merge TPM values
pd2 <- tpm.long[goi]
pd2[, gene.symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
pd2[is.na(gene.symbol), gene.symbol := gene]
expr <- pd2[, .(any(call == TRUE)), by = gene][V1 == TRUE, unique(gene)] # don't do this, use list of expressed genes instead!!!
pd2 <- pd2[expr]

# weird plot
ggplot(pd2,
       aes(x = accession, y = tpm, shape = accession, colour = stage,
           group = accession)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30),
        axis.ticks.length = unit(0, "mm")) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ gene.symbol, scales = "free_y") +
  geom_smooth(mapping = aes(group = stage), method = "loess", se = FALSE) +
  geom_point(position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1))
