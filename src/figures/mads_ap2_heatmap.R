#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(cowplot)
library(data.table)
library(DESeq2)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)

###########
# GLOBALS #
###########

tfdb_file <- snakemake@input[["tfdb"]]
families_file <- snakemake@input[["families"]]
vst_file <- snakemake@input[["vst"]]
leading_edge_file <- snakemake@input[["leading_edge"]]
arora_file <- snakemake@input[["arora"]]
arora_subclades_file <- snakemake@input[["arora_subclades"]]
sharoni_file <- snakemake@input[["sharoni"]]

# plots
fig1_file <- snakemake@output[["fig1"]]
sf1_file <- snakemake@output[["sf1"]]


# dev
# tfdb_file <- "output/010_data/tfdb.Rds"
# families_file <- "output/010_data/tfdb_families.Rds"
# vst_file <- "output/050_deseq/vst.Rds"
# arora_file <- "data/genome/os/arora.csv"
# arora_subclades_file <- "data/genome/os/arora_subclades.csv"
# sharoni_file <- "data/genome/os/sharoni_table_s1.csv"
# leading_edge_file <- "tmp/no_pca/leading_edge.csv"
# wald_stage_file <- "output/050_deseq/wald_tests/expr_genes/all/stage.csv"
# fig1_file <- "tmp/no_pca/mads_ap2_heatmap.pdf"

gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

spec_order <- c("or" = "O. rufipogon",
                "osi" = "O. sativa indica",
                "ob" = "O. barthii",
                "og" = "O. glaberrima")

stage_order <- c(PBM = "IM", SM = "DM")

########
# MAIN #
########

# data
vst <- readRDS(vst_file)
tfdb <- readRDS(tfdb_file)
families <- readRDS(families_file)
arora <- fread(arora_file)
arora_subclades <- fread(arora_subclades_file)
sharoni <- fread(sharoni_file)
leading_edge <- fread(leading_edge_file)
wald_stage <- fread(wald_stage_file)

#################
# Heatmap panel #  
#################

# protect the class column, if it's not a MADS gene
families[!is.na(Class) & Family != "MADS",
         Class := paste0('"', Class, '"')]

# use Arora classes for MADS genes
families[Family == "MADS",
         Class := arora[TIGR == `Protein ID`, unique(type)],
         by = `Protein ID`]
families[Family == "MADS" & grepl("_", Class), Class := NA]

# add subclades where possible
arora_subclades[!is.na(arora_subclade) & !grepl("-like", arora_subclade),
                arora_subclade := paste0('italic("', arora_subclade, '")')]
arora_subclades[!is.na(arora_subclade),
                arora_subclade := gsub("([[:alnum:]]+)-like",
                                       'italic("\\1-")*"like"',
                                       arora_subclade)]
arora_subclade_genes <- arora_subclades[!is.na(arora_subclade),
                                        unique(gene_id)]
families[`Protein ID` %in% arora_subclade_genes,
         Class := arora_subclades[gene_id == `Protein ID`,
                                  unique(arora_subclade)],
         by = `Protein ID`]

# use sharoni clades for AP2. it's a nasty mess
sharoni_classes <- sharoni[!is.na(`MSU locus ID`) &
                             `MSU locus ID` != "" &
                             grepl("^Os", `MSU locus ID`),
                           .(gene_id = paste("LOC", `MSU locus ID`, sep = "_"),
                             Class = `Phy. Subfamily`)]
sharoni_genes <- sharoni_classes[, unique(gene_id)]
families[`Protein ID` %in% sharoni_genes,
         Class := paste0('"',
                         sharoni_classes[gene_id == `Protein ID`, unique(Class)],
                         '"'),
         by = `Protein ID`]

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

# map stage
mean_vst[, stage := factor(plyr::revalue(stage, stage_order),
                           levels = stage_order)]

# scale & centre transformed counts
mean_vst[, scaled_vst := scale(mean_vst), by = gene_id]

# select genes to plot
wald_stage[, alfc := abs(log2FoldChange)]
setkey(wald_stage)
setorder(wald_stage, -alfc, na.last = TRUE)
wald_stage[, lfc_rank := order(alfc, decreasing = TRUE, na.last = TRUE)]
top10pct <- wald_stage[lfc_rank / max(lfc_rank) < 0.1,
           unique(gene_id)]
plot_mads <- unique(intersect(mads_genes, top10pct))
plot_ap2 <- unique(intersect(ap2_genes, top10pct))

# lfc_cutoff <- 1
# leading_edge <- unique(leading_edge, by = "gene_id")
# plot_ap2 <- leading_edge[abs(log2FoldChange) > lfc_cutoff &
#                            family == "AP2-EREBP", unique(gene_id)]
# plot_mads <- leading_edge[abs(log2FoldChange) > lfc_cutoff &
#                             family == "MADS", unique(gene_id)]

# set up scale
v_max <- mean_vst[gene_id %in% c(plot_ap2, plot_mads), max(abs(scaled_vst))]
v_lim <- c(-v_max,
           v_max)

# for now this relies on the objects in the environment
PlotHeatmapWithFamily <- function(plot_genes, plot_title) {
  # plot_genes <- plot_ap2
  # plot_title <- "AP2"
  
  # cut by lfc sign
  pd <- mean_vst[gene_id %in% plot_genes & species %in% names(spec_order)]
  pd <- merge(pd,
              wald_stage[, .(gene_id, log2FoldChange)],
              by = "gene_id")
  pd[, cut_row := log2FoldChange < 0]
  
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
    xlab(" ") + # dummy xlab prevents axis text being cut off in some versions
    ylab(NULL) + 
    ggtitle(plot_title) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_viridis_c(
      #limits = v_lim,
      guide = guide_colourbar(title = "Scaled\nreads")) +
    facet_grid(cut_row ~ stage, scales = "free_y", space = "free_y") +
    geom_raster()
  
  # check number of classes
  pd[, family := families[`Protein ID` == gene_id, Class[[1]]], by = gene_id]
  if(pd[, length(unique(family))] == 1){
    return(hm)
  }
  
  # for multi-class families, add the class column
  pd[, family := factor(family,
                        levels = sort(unique(as.character(family))))]
  pd[is.na(family), family := '"Other"'] 
  famplot <- ggplot(pd, aes(y = label, x = "Clade", fill = family)) +
    theme_minimal(base_size = 8, base_family = "Helvetica") +
    theme(strip.text.y = element_blank(),
          axis.text = element_blank(),
          axis.text.x = element_blank(),
          legend.key.size = unit(0.8, "lines"),
          panel.grid = element_blank()) +
    xlab(NULL) + ylab(NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_brewer(palette = "Set3",
                      guide = guide_legend(title = "Clade",
                                           label.hjust = 0,
                                           label.vjust = 0.5),
                      labels = function(l) parse(text = l)) +
    facet_grid(cut_row ~ "Clade", scales = "free_y", space = "free_y") +
    geom_raster()
  
  # combine the plots
  hmg <- ggplotGrob(hm)
  famplotg <- ggplotGrob(famplot)
  
  # get the famplot bits
  famplot_panel <- gtable_filter(famplotg, "panel")
  famplot_strip <- gtable_filter(famplotg, "strip-t-1")
  
  # make room
  hmg2 <- gtable_add_cols(hmg,
                          unit(c(1.25/4, 4), c("null", "pt")),
                          4)
  
  # insert it into the table
  hmg3 <- gtable_add_grob(hmg2, famplot_panel, 8, 5, 10, 5 )
  hmg4 <- gtable_add_grob(hmg3, famplot_strip, 7, 5, 7, 5)
  
  # get the legends
  famplot_legend <- gtable_filter(famplotg, "guide-box")
  hmg_legend <- gtable_filter(hmg, "guide-box")
  
  # unify the legend widths
  hmg_legend$widths <- famplot_legend$widths
  
  # combine the legends
  both_legends <- gtable_matrix(name = "legends",
                                grobs = matrix(list(hmg_legend, famplot_legend),
                                               ncol = 1),
                                widths = hmg_legend$widths,
                                heights = unit(c(1, 1), "null"))
  
  # add room between legends
  both_legends$heights <- unit(c(0.3, 0.3), "npc") # specify heights
  both_legends2 <- gtable_add_row_space(both_legends, unit(0.05, "npc")) # between
  both_legends3 <- gtable_add_rows(both_legends2, unit(0.15, "npc"), 0) # top
  both_legends4 <- gtable_add_rows(both_legends3, unit(0.2, "npc"), -1) # bottom
  
  # put the legend in the main plot
  hmg5 <- gtable_remove_grobs(hmg4, "guide-box")
  hmg6 <- gtable_add_grob(hmg5, both_legends4, 8, 14, 10, 14,
                          name = "both_legends")
  
  # add width for the legends
  hmg6$widths[[14]] <- famplot_legend$widths
  
  return(hmg6)
}

ap2_gt <- PlotHeatmapWithFamily(plot_ap2,
                                expression(italic("AP2/EREBP-")*"like"))
mads_gt <- PlotHeatmapWithFamily(plot_mads,
                                 expression(italic("MADS-")*"box"))

###########
# COMBINE #
###########

cowplot <- plot_grid(ap2_gt,
                     mads_gt,
                     ncol = 2,
                     label_size = 10,
                     label_fontfamily = "Helvetica",
                     rel_widths = c(0.9, 1))

ggsave(fig1_file,
       device = cairo_pdf,
       cowplot,
       width = 178,
       height = 150,
       units = "mm")


#################################
# SUPPLEMENTARY HOMEOBOX FIGURE #
#################################

plot_hb <- leading_edge[abs(log2FoldChange) > lfc_cutoff &
                          family == "HB", unique(gene_id)]
hb_gt <- PlotHeatmapWithFamily(plot_hb, "Homeobox")

plot_nac <- leading_edge[abs(log2FoldChange) > lfc_cutoff &
                           family == "NAC", unique(gene_id)]
nac_gt <- PlotHeatmapWithFamily(plot_nac, "NAC")

plot_spl <- leading_edge[abs(log2FoldChange) > lfc_cutoff &
                           family == "SBP", unique(gene_id)]
spl_gt <- PlotHeatmapWithFamily(plot_spl, "SBP")

# combine it
s8_panels <- plot_grid(
  hb_gt, nac_gt, spl_gt,
  ncol = 3,
  rel_widths = c(1.1, 1, 1))

ggsave(sf1_file,
       device = cairo_pdf,
       s8_panels,
       width = 178*3/2,
       height = 150,
       units = "mm")

# Log
sessionInfo()

