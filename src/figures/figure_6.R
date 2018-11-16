#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(cowplot)

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
    theme_minimal(base_size = 8, base_family = "Helvetica") +
    theme(axis.text.x = element_text(angle= 30,
                                     face = "italic",
                                     hjust = 1,
                                     vjust = 1),
          strip.text = element_text(face = "italic"),
          panel.border = element_rect(fill = NA,
                                      colour = "black")) +
    scale_fill_manual(values = bar_fills, guide = FALSE) +  
    xlab(NULL) +
    ylab(expression("L"[2]*"FC" %+-% "SE")) +
    facet_wrap(~label, scales = "free_y", ncol = ncol) +
    geom_col(width = 0.8) +
    geom_hline(yintercept = 0, size = 0.2) +
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

acc_colours <- RColorBrewer::brewer.pal(4, "Paired")[c(1, 2, 3, 4)]
wald_file <- "output/050_deseq/wald_tests/tfs/all/stage_within_species.csv"
domestication_file <- "output/050_deseq/wald_tests/tfs/all/domestication.csv"
indica_file <- "output/050_deseq/wald_tests/tfs/all/stage_accession_indica.csv"
glab_file <- "output/050_deseq/wald_tests/tfs/all/stage_accession_glaberrima.csv"

########
# MAIN #
########

stage_lfc <- fread(wald_file)
domestication <- fread(domestication_file)
indica <- fread(indica_file)
glab <- fread(glab_file)

# mung the stage_lfc
stage_lfc[, accession := gsub(".*group ([[:alpha:]]+).*", "\\1", wald_test)]
stage_lfc[, accession := factor(plyr::revalue(accession, spec_order),
                                levels = spec_order)]

# don't plot japonica
plot_lfc <- stage_lfc[accession != "O. sativa japonica"]

# lists of genes
dom_genes <- domestication[padj < alpha, unique(gene_id)]
ind_genes <- indica[padj < alpha, unique(gene_id)]
glab_genes <- glab[padj < alpha, unique(gene_id)]
both_sep <- intersect(ind_genes, glab_genes)

# make the plots
dom_gp <- PlotGoiLfcs(dom_genes,
                      lfc_table = plot_lfc,
                      ncol = 5)

sep_but_both <- PlotGoiLfcs(
  both_sep[!both_sep %in% dom_genes],
  lfc_table = plot_lfc,
  ncol = 6)


cowplot <- plot_grid(dom_gp, sep_but_both,
          ncol = 1,
          align = "hv",
          axis = "tblr",
          labels = c("A", "B"),
          label_size = 10,
          label_fontfamily = "Helvetica",
          rel_heights = c(4, 1.25))

ggsave("test/Figure_6.pdf",
       cowplot,
       width = 178,
       height = 178,
       units = "mm")

