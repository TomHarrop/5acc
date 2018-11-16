library(cowplot)
library(data.table)
library(DESeq2)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)

gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

vst_file <- "output/050_deseq/vst.Rds"
tfdb_file <- "output/010_data/tfdb.Rds"
# test with domestication genes
dom_genes_file <- "output/050_deseq/wald_tests/tfs/sig/domestication.csv"
pcro_file <- "output/050_deseq/rlog_pca/pcro.Rds"

spec_order <- c("or" = "O. rufipogon",
                "osi" = "O. sativa indica",
                "osj" = "O. sativa japonica",
                "ob" = "O. barthii",
                "og" = "O. glaberrima")

########
# MAIN #
########

# data
vst <- readRDS(vst_file)
dom_genes <- fread(dom_genes_file)
pcro <- data.table(readRDS(pcro_file))
tfdb <- readRDS(tfdb_file)

# list of genes
ap2_genes <- tfdb[Family == "AP2-EREBP", unique(`Protein ID`)]
mads_genes <- tfdb[Family == "MADS", unique(`Protein ID`)]

# mean vst
vst_wide <- as.data.table(assay(vst), keep.rownames = TRUE)
setnames(vst_wide, "rn", "gene_id")
vst_dt <- melt(vst_wide,
               id.vars = "gene_id",
               value.name = "vst",
               variable.name = "library")
vst_dt[, c("species", "stage", "rep") := tstrsplit(library, "_")]
mean_vst <- vst_dt[, .(mean_vst = gm_mean(vst)),
                   by = .(gene_id, species, stage)]

# scale & centre transformed counts
mean_vst[, scaled_vst := scale(mean_vst), by = gene_id]

# select genes to plot
pcro[, abs_pc5 := abs(PC5)]
setorder(pcro, abs_pc5, na.last = TRUE)
pcro[!is.na(abs_pc5), abs_pc5_rank := order(abs_pc5, decreasing = TRUE)]
top_10pct <- pcro[abs_pc5_rank / max(abs_pc5_rank) < 0.1]
plot_ap2 <- intersect(ap2_genes,
                      top_10pct[, unique(locus_id)])
plot_mads <- intersect(mads_genes,
                       top_10pct[, unique(locus_id)])

# set up scale
v_max <- mean_vst[gene_id %in% c(plot_ap2, plot_mads), max(abs(scaled_vst))]
v_lim <- c(-v_max,
           v_max)


# for now this relies on the objects in the environment
PlotHeatmapWithFamily <- function(plot_genes, plot_title) {
  # plot_genes <- plot_ap2
  # plot_title <- "AP2"
  
  # cut by posn on PC5
  pd <- mean_vst[gene_id %in% plot_genes]
  pd <- merge(pd,
              pcro[, .(locus_id, PC5)],
              by.x = "gene_id", by.y = "locus_id")
  pd[, cut_row := PC5 < 0]
  
  # label genes
  pd[, symbol := oryzr::LocToGeneName(gene_id)$symbols[[1]], by = gene_id]
  pd[, label := ifelse(is.na(symbol),
                       gene_id,
                       paste(symbol, gene_id, sep = "\n")),
     by = gene_id]
  
  # order genes
  hc_m <- as.matrix(data.frame(
    dcast(pd,
          label ~ species + stage,
          value.var = "mean_vst"),
    row.names = "label"))
  hc <- hclust(dist(hc_m, method = "minkowski"))
  gene_order <- hc$labels[hc$order]
  pd[, label := factor(label, levels = gene_order)]
  
  # order species
  pd[, species := factor(plyr::revalue(species, spec_order),
                         levels = spec_order)]
  
  # ggplot heatmap
  hm <- ggplot(pd, aes(y = label, x = species, fill = scaled_vst)) +
    theme_minimal(base_size = 8, base_family = "Helvetica") +
    theme(strip.text.y = element_blank(),
          axis.text = element_text(face = "italic"),
          axis.text.y = element_text(size = 4),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          plot.title = element_text(hjust = 0.5),
          legend.key.size = unit(0.8, "lines"),
          legend.justification = "left") +
    xlab(NULL) + ylab(NULL) + ggtitle(plot_title) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_viridis_c(
      #limits = v_lim,
      guide = guide_colourbar(title = "Scaled\nreads")) +
    facet_grid(cut_row ~ stage, scales = "free_y", space = "free_y") +
    geom_raster()
  
  # family panel, dummy for now
  pd[, family := sample(letters[1:2])[[1]], by = gene_id]
  famplot <- ggplot(pd, aes(y = label, x = "Type", fill = family)) +
    theme_minimal(base_size = 8, base_family = "Helvetica") +
    theme(strip.text.y = element_blank(),
          axis.text = element_blank(),
          axis.text.x = element_blank(),
          legend.key.size = unit(0.8, "lines"),
          legend.justification = "left") +
    xlab(NULL) + ylab(NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_brewer(palette = "Set1",
                      guide = guide_legend(title = "Type")) +
    facet_grid(cut_row ~ "Type", scales = "free_y", space = "free_y") +
    geom_raster()

  # combine the plots
  hmg <- ggplotGrob(hm)
  famplotg <- ggplotGrob(famplot)
  
  # get the famplot bits
  famplot_panel <- gtable_filter(famplotg, "panel")
  famplot_strip <- gtable_filter(famplotg, "strip-t-1")
  
  # make room
  hmg2 <- gtable_add_cols(hmg, unit(c(1/5, 4), c("null", "pt")),4)
  
  # insert it into the table
  hmg3 <- gtable_add_grob(hmg2, famplot_panel, 8, 5, 10, 5 )
  hmg4 <- gtable_add_grob(hmg3, famplot_strip, 7, 5, 7, 5)
  
  # get the legends
  famplot_legend <- gtable_filter(famplotg, "guide-box")
  hmg_legend <- gtable_filter(hmg, "guide-box")
  
  famplot_legend$widths <- hmg_legend$widths
  
  both_legends <- gtable_matrix(name = "legends",
                                grobs = matrix(list(hmg_legend, famplot_legend),
                                               ncol = 1),
                                widths = hmg_legend$widths,
                                heights = unit(c(1, 1), "null"))
  both_legends2 <- gtable_add_rows(both_legends, unit(0.5, "null"), 1)
  both_legends3 <- gtable_add_rows(both_legends2, unit(2, "null"), -1)
  both_legends4 <- gtable_add_rows(both_legends3, unit(2, "null"), 0)
  
  hmg5 <- gtable_remove_grobs(hmg4, "guide-box")
  hmg6 <- gtable_add_grob(hmg5, both_legends4, 8, 14, 10, 14)
  
  return(hmg6)
}

ap2_gt <- PlotHeatmapWithFamily(plot_ap2, expression(italic("AP2/EREBP")*"-like"))
mads_gt <- PlotHeatmapWithFamily(plot_mads, expression(italic("MADS")))

cowplot <- plot_grid(ap2_gt,
                     mads_gt,
                     ncol = 2,
                     labels = "C",
                     label_size = 10,
                     label_fontfamily = "Helvetica")
ggsave("test/Figure_3C.pdf",
       device = cairo_pdf,
       cowplot,
       width = 178,
       height = 150,
       units = "mm")
