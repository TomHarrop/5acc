#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)
library(ggtree)
library(grid)
library(gridExtra)
library(cowplot)


###########
# GLOBALS #
###########

tpm_file <- snakemake@input[["tpm"]]
cluster_file <- snakemake@input[["cluster"]]
wald_file <- snakemake@input[["wald"]]
hyperg_file <- snakemake@input[["hyperg"]]
correlation_file <- snakemake@input[["correlation"]]
pheno_key_file <- snakemake@input[["pheno_key"]]
cali_corr_file <- snakemake@input[["cali_corr"]]

# plots
fig1_file <- snakemake@output[["fig1"]]
sf1_file <- snakemake@output[["sf1"]]


# dev
# tpm_file <- "output/060_tpm/tpm_with_calls.Rds"
# cluster_file <- "output/070_clustering/tfs/annotated_clusters_scaled_l2fc.csv"
# wald_file <- "output/050_deseq/wald_tests/tfs/all/stage_within_species.csv"
# hyperg_file <- "output/070_clustering/tfs/hypergeom.csv"
# correlation_file <- "output/080_phenotype/mtp_cluster_correlation.csv"
# pheno_key_file <- "data/phenotyping/phenotype_name_key.csv"
# cali_corr_file <- "output/080_phenotype/cali_cluster_correlation.csv"


spec_order <- c("rufipogon" = "O. rufipogon",
                "indica" = "O. sativa indica",
                "barthii" = "O. barthii",
                "glaberrima" = "O. glaberrima")

# plot colours
base_colour = "#FFFFFF"
rdbu <- rev(RColorBrewer::brewer.pal(5, "RdBu"))
prgn <- RColorBrewer::brewer.pal(5, "PRGn")
ylgnbu <- RColorBrewer::brewer.pal(5, "YlGnBu")
Purples <- RColorBrewer::brewer.pal(5, "Purples")
set1 <- RColorBrewer::brewer.pal(9, "Set1")

########
# MAIN #
########

# read data
clusters <- fread(cluster_file)
tpm <- readRDS(tpm_file)
wald_results <- fread(wald_file)
hyperg_results <- fread(hyperg_file)
correlations <- fread(correlation_file)
pheno_key <- fread(pheno_key_file)
cali_corr <- fread(cali_corr_file)

# set up a theme for the heatmaps
theme_hm <- theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 8),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(colour = "black"),
        plot.margin = unit(c(1, 2, 1, 0), "mm"),
        legend.position = "top")

# average expression value per cluster, should be equivalent to cluster core
clusters_long <- melt(clusters,
                      id.vars = c("MsuID", "cluster"),
                      measure.vars = names(spec_order),
                      variable.name = "accession",
                      value.name = "scaled_expression_value")

cluster_heatmap <- clusters_long[
  , .(core_value = mean(scaled_expression_value)),
  by = .(cluster, accession)]

# generate a tree for the cluster cores
core_matrix <- as.matrix(
  data.frame(
    dcast(cluster_heatmap,
          cluster ~ accession,
          value.var = "core_value"),
    row.names = "cluster"))
hc <- hclust(dist(core_matrix,
                  method = "minkowski"),
             method = "ward.D")
phylo <- ape::as.phylo(hc)

# generate a heatmap for cores
cluster_order <- hc$order
cluster_heatmap[, cluster := factor(cluster, levels = cluster_order)]
cluster_heatmap[, accession := factor(plyr::revalue(accession, spec_order),
                                      levels = spec_order)]

gp <- ggplot(cluster_heatmap, aes(y = cluster, x = accession, fill = core_value)) +
  theme_hm +
  theme(axis.text.x = element_text(face = "italic"),
        legend.position = "right",
        legend.key.size = unit(0.8, "lines")) +
  scale_y_discrete(breaks = c(1:7), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colours = rdbu,
                       guide = guide_colourbar(title = NULL)) +
                         #title = expression("L"[2] * "FC"),
                         #title.position = "top",
                         #title.hjust = 0.5)) +
  geom_raster() +
  ggtitle(expression("L"[2]*"FC"))
gp

# generate a heatmap for the hypergeom
families <- hyperg_results[p_adj < 0.1, unique(family)]
hyperg_to_plot <- hyperg_results[family %in% families]
hyperg_to_plot[, enrichment := (n_in_cluster/all_in_cluster)/(n_in_bg / all_in_bg)]
hyperg_to_plot[, cluster := factor(cluster, levels = cluster_order)]

# format p-values for printint
hyperg_to_plot[, label := format(p_adj,
                                 digits = 2,
                                 nsmall = 1,
                                 scientific = TRUE,
                                 drop0trailing = FALSE)]
hyperg_to_plot[!is.na(p_adj),
               label := paste0(
                 # "italic(p)['adj'] ==",
                 gsub("^([^e]+).*", "'\\1'", label),
                 " %*% 10^",
                 gsub(".*e([^e]+)$", "\\1", label))]
hyperg_to_plot[p_adj < 0.1, label2 := "'*'"]

enrichment_plot <- ggplot(hyperg_to_plot, aes(x = family,
                                              y = cluster,
                                              fill = enrichment,
                                              label = label2)) +
  theme_hm +
  scale_y_discrete(breaks = c(1:7), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colours = ylgnbu,
                       na.value = base_colour,
                       guide = guide_colourbar(title = "TF enrichment",
                                               title.position = "top",
                                               title.hjust = 0.5)) +
  geom_raster() +
  geom_text(parse = TRUE)
enrichment_plot

# phenotype correlations
corr_vars <- c("primary_branch_number", "secondary_branch_number", "spikelet_number")
correlations[, cluster := factor(cluster, levels = cluster_order)]
correlations[, variable := factor(variable, levels = corr_vars)]
correlations_pd <- merge(correlations,
                         pheno_key,
                         by.x = "variable",
                         by.y = "full_name",
                         all.x = TRUE,
                         all.y = FALSE)

corr_plot <- ggplot(correlations_pd, aes(x = short_name,
                                         y = cluster,
                                         fill = pearson_correlation)) +
  theme_hm +
  scale_y_discrete(breaks = c(1:7), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colours = rdbu,
                       guide = guide_colourbar(title = "Pearson correlation",
                                               title.position = "top",
                                               title.hjust = 0.5)) +
  geom_raster()
corr_plot

# horizontal bar plot
spn_corr <- ggplot(correlations_pd[short_name == "SpN"],
                    aes(y = pearson_correlation, x = cluster, fill = pearson_correlation)) +
  theme_hm +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) +
  coord_flip() +
  scale_x_discrete(breaks = c(1:7), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(guide = FALSE) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  ylim(c(-1, 1)) +
  ggtitle("SpN")
spn_corr

sbn_corr <- ggplot(correlations_pd[short_name == "SBN"],
                     aes(y = pearson_correlation, x = cluster, fill = pearson_correlation)) +
  theme_hm +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) +
  coord_flip() +
  scale_x_discrete(breaks = c(1:7), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(guide = FALSE) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  ylim(c(-1, 1)) +
  ggtitle("SBN")




pbn_corr <- ggplot(correlations_pd[short_name == "PBN"],
                     aes(y = pearson_correlation, x = cluster, fill = pearson_correlation)) +
  theme_hm +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) +
  coord_flip() +
  scale_x_discrete(breaks = c(1:7), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(guide = guide_colourbar(title = "Correlation with PBN",
                                               title.position = "top",
                                               title.hjust = 0.5)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  ylim(c(-1, 1))

cali_corr[, cluster := factor(cluster, levels = cluster_order)]
pc1_corr <- ggplot(cali_corr[PC == "PC1"],
                     aes(y = pearson_correlation, x = cluster, fill = pearson_correlation)) +
  theme_hm +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) +
  coord_flip() +
  scale_x_discrete(breaks = c(1:7), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(guide = FALSE) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  ylim(c(-1, 1)) +
  ggtitle("PC1")
pc1_corr

# number of genes
gene_no <- ggplot(clusters[, length(unique(MsuID)), by = cluster],
                  aes(y = V1, x = cluster, fill = V1)) +
  theme_hm +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) +
  coord_flip() +
  scale_x_discrete(breaks = c(1:7), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn(colours = Purples,
                       guide = guide_colourbar(title = "Genes",
                                               title.position = "top",
                                               title.hjust = 0.5)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0)
gene_no

# generate a tree
p <- ggtree(phylo) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_blank(),
        plot.margin = unit(c(0, 0, 0, 2), "mm")) +
  ylab("Cluster") +
  scale_y_discrete(breaks = c(1:7)) +
  scale_x_continuous(expand = c(0, 0.5)) +
  geom_tiplab()
p

# join plots
cowplot <- plot_grid(p,
                     gp,
                     pc1_corr,
                     spn_corr,
                     sbn_corr,
                     nrow = 1,
                     align = "h",
                     axis = "tb",
                     rel_widths = c(1.3, 4, 2, 2, 2),
                     labels = c("A", "", "B", "", ""),
                     label_size = 10,
                     label_fontfamily = "Helvetica",
                     label_x = c(0, -0.1, -0.1, 0, 0))
ggsave(fig1_file,
       device = cairo_pdf,
       cowplot,
       width = 178,
       height = 75,
       units = "mm")

############################
# PLOT INDIVIDUAL CLUSTERS #
############################

# extract per-accession stats
wald_results[, accession := gsub(".*group ([[:alpha:]]+).*", "\\1", wald_test)]
acc_results <- wald_results[, .(accession, gene_id, log2FoldChange, lfcSE, symbols)]
acc_results[!symbols == "", label := paste(gene_id, symbols, sep = " | ")]
acc_results[symbols == "", label := gene_id]

# plot cluster 5 only
plot_genes <- clusters[cluster == 5, unique(MsuID)]
my_fills <- RColorBrewer::brewer.pal(4, "Paired")
pd2 <- acc_results[gene_id %in% plot_genes & accession %in% names(spec_order)]
pd2[, accession := factor(plyr::revalue(accession, spec_order),
                          levels = spec_order)]
gp <- ggplot(pd2, aes(x = accession,
                      y = log2FoldChange,
                      ymin = log2FoldChange - lfcSE,
                      ymax = log2FoldChange + lfcSE,
                      fill = accession)) +
  theme_grey(base_size = 6) +
  theme(axis.text.x =  element_text(face = "italic"),
        strip.text = element_text(face = "italic")) +
  xlab(NULL) +
  ylab(expression(L[2] * "FC ± SE (IM vs. DM)")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_manual(values = my_fills, guide = FALSE) +  
  facet_wrap(~ label, scales = "free_y", nrow = 4) +
  geom_col(position = "dodge", alpha = 0.5) +
  geom_errorbar(width = 0.1) 


ggsave(sf1_file,
       gp,
       device = cairo_pdf,
       width = 178,
       height = 225*2/3,
       units = "mm")


# Log
sessionInfo()

