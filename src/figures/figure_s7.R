#!/usr/bin/env Rscript

# Figure S7: Expression of *AP2/EREBP*-like genes in *O. sativa japonica* cv.
# Nipponbare meristems [data from @harropGeneExpressionProfiling2016].

library(data.table)
library(ggplot2)

tpm_data_file <- "data/tpj13147-sup-0007-datas1.tsv"
goi <- c("LOC_Os07g03250", "LOC_Os05g32270")
libs <- c("n1r1",
          "n1r3",
          "n1r4",
          "n2r1",
          "n2r3",
          "n2r4",
          "n3r1",
          "n3r2", 
          "n3r3",
          "n4r1",
          "n4r2",
          "n4r3")

# extract plot data
tpm_data <- as.data.table(readr::read_tsv(tpm_data_file))
goi_data <- tpm_data[msuId %in% goi]
goi_pd <- melt(goi_data,
     id.vars = "msuId",
     measure.vars = libs,
     variable.name = "library",
     value.name = "Expression (TPM)")

# add stage
goi_pd[, stage := substr(library, start = 1, stop = 2)]
old <- c("n1", "n2", "n3", "n4")
new <- c("RM", "PBM", "ePBM/AM", "SM")
goi_pd[, stage := factor(plyr::mapvalues(stage, from = old, to = new),
                         levels = new)]

# add gene labels
goi_pd[, symbol := oryzr::LocToGeneName(msuId)$symbols]

# set up labels and order by msuId
goi_pd[!is.na(symbol), symbol := paste(msuId, "|", symbol)]
goi_pd[is.na(symbol), symbol := msuId]
goi_pd[, symbol := factor(symbol, levels = unique(symbol))]

# plot
gp <- ggplot(goi_pd, aes(x = stage, y = `Expression (TPM)`, group = symbol)) +
    theme_minimal(base_size = 8, base_family = "Helvetica") +
    theme(strip.text = element_text(face = "italic"),
          panel.background = element_rect(colour = "black")) +
    facet_wrap(~symbol) +
    xlab(NULL) +
    ylab("Expression (TPM)") +
    stat_smooth(se = FALSE, colour = "grey", size = 0.5) + 
    geom_point(shape = 16, position = position_jitter(width = 0.2))

ggsave("figures/Figure_S7.pdf",
       gp,
       device = cairo_pdf,
       width = 114,
       height = 50,
       units = "mm")

