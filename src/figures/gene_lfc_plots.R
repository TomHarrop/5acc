#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

#############
# FUNCTIONS #
#############

PlotGoiLfcs <- function(genes_to_plot,
                        lfc_table = stage_lfc,
                        bar_fills = acc_colours,
                        ncol = 6) {
  gene_pd <- lfc_table[gene_id %in% genes_to_plot]
  gene_pd[!symbols == "", label := paste(gene_id, symbols, sep = "\n")]
  gene_pd[symbols == "", label := gene_id]
  gp <- ggplot(gene_pd, aes(x = accession,
                            y = log2FoldChange,
                            ymin = log2FoldChange - lfcSE,
                            ymax = log2FoldChange + lfcSE,
                            fill = accession)) +
    theme_grey(base_size = 6, base_family = "Helvetica") +
    theme(axis.text.x = element_text(angle= 30,
                                     face = "italic",
                                     hjust = 1,
                                     vjust = 1),
          strip.text = element_text(face = "italic")) +
    scale_fill_manual(values = bar_fills, guide = FALSE) +  
    xlab(NULL) +
    ylab(expression("L"[2]*"FC" %+-% "SE")) +
    facet_wrap(~label, scales = "free_y", ncol = ncol) +
    geom_col(width = 0.8) +
    geom_errorbar(width = 0.2, alpha = 0.5)
  return(gp)
}

########
# DATA # 
########

alpha <- 0.1

spec_order <- c("rufipogon" = "O. rufipogon",
                "indica" = "O. sativa indica",
                "japonica" = "O. sativa japonica",
                "barthii" = "O. barthii",
                "glaberrima" = "O. glaberrima")

acc_colours <- RColorBrewer::brewer.pal(4, "Paired")[c(1, 2, 2, 3, 4)]

tpm_file <- "output/060_tpm/tpm_with_calls.Rds"
wald_file <- "output/050_deseq/wald_tests/tfs/all/stage_within_species.csv"
domestication_file <- "output/050_deseq/wald_tests/tfs/all/domestication.csv"
indica_file <- "output/050_deseq/wald_tests/tfs/all/stage_accession_indica.csv"
glab_file <- "output/050_deseq/wald_tests/tfs/all/stage_accession_glaberrima.csv"

########
# MAIN #
########

tpm <- readRDS(tpm_file)
stage_lfc <- fread(wald_file)
domestication <- fread(domestication_file)
indica <- fread(indica_file)
glab <- fread(glab_file)

# mung the stage_lfc
stage_lfc[, accession := gsub(".*group ([[:alpha:]]+).*", "\\1", wald_test)]
stage_lfc[, accession := factor(plyr::revalue(accession, spec_order),
                                levels = spec_order)]

# make the plots
dom_gp <- PlotGoiLfcs(domestication[padj < alpha, unique(gene_id)],
                      ncol = 4)
ggsave("test/dom_both.pdf",
       dom_gp, 
       width = 190,
       height = 277,
       units = "mm")

indica_gp <- PlotGoiLfcs(indica[padj < alpha, unique(gene_id)],
                         ncol = 8)
ggsave("test/indica_dom.pdf",
       indica_gp, 
       width = 190,
       height = 277,
       units = "mm")

glab_gp <- PlotGoiLfcs(glab[padj < alpha, unique(gene_id)],
                       ncol = 6)
ggsave("test/glab_dom.pdf",
       glab_gp, 
       width = 190,
       height = 277,
       units = "mm")

sep_but_both <- PlotGoiLfcs(
  intersect(indica[padj < alpha, unique(gene_id)],
            glab[padj < alpha, unique(gene_id)]),
  ncol = 4)
ggsave("test/sep_but_both.pdf",
       sep_but_both, 
       width = 190,
       height = 277,
       units = "mm")

