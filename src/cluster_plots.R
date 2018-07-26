library(data.table)
library(ggplot2)
library(ggtree)
library(grid)
library(gridExtra)

tpm_file <- "output/060_tpm/tpm_with_calls.Rds"
cluster_file <- "output/070_clustering/tfs/annotated_clusters_scaled_l2fc.csv"
wald_file <- "output/050_deseq/wald_tests/tfs/all/stage_within_species.csv"
hyperg <- "output/070_clustering/tfs/hypergeom.csv"

#############
# FUNCTIONS #
#############

MakeFakeLabels <- function(x){
  rep(accession_order[[which.max(nchar(accession_order))]], length(x))
}

# what about ALOGs (not annotated on TFDB)
ALOG <- c('LOC_Os07g04670', 'LOC_Os02g07030', 'LOC_Os06g46030',
          'LOC_Os02g41460', 'LOC_Os04g43580', 'LOC_Os10g33780',
          'LOC_Os02g56610', 'LOC_Os01g61310', 'LOC_Os05g39500',
          'LOC_Os05g28040')

###########
# GLOBALS #
###########

accession_order <- c("rufipogon", "indica", "barthii",  "glaberrima")

# plot colours
base_colour = "#FFFFFF"
rdbu <- RColorBrewer::brewer.pal(7, "RdBu")
set1 <- RColorBrewer::brewer.pal(9, "Set1")

########
# MAIN #
########

# read data
clusters <- fread(cluster_file)
tpm <- readRDS(tpm_file)
wald_results <- fread(wald_file)
hyperg_results <- fread(hyperg)

# average expression value per cluster, should be equivalent to cluster core
clusters_long <- melt(clusters,
                      id.vars = c("MsuID", "cluster"),
                      measure.vars = accession_order,
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
cluster_heatmap[, accession := factor(accession, levels = accession_order)]

gp <- ggplot(cluster_heatmap, aes(y = cluster, x = accession, fill = core_value)) +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        panel.background = element_rect(colour = "black")) +
  ggtitle("Core expression") +
  scale_y_discrete(breaks = c(1:7), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(low = set1[[2]],
                       mid = base_colour, 
                       high = set1[[1]],
                       guide = guide_colourbar(title = expression("Scaled" ~ "L"[2] * "FC"))) +
  geom_raster()
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
  theme_void(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        panel.background = element_rect(colour = "black")) +
  ggtitle("TF enrichment") +
  scale_y_discrete(breaks = c(1:7), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colours = c(base_colour, set1[[4]]),
                       guide = guide_colourbar(title = "Enrichment")) +
  geom_raster() +
  geom_text(parse = TRUE)
enrichment_plot


# phenotype correlations
dummy_plot <- ggplot(data.frame(x = 1, y = 1), aes(x,y)) +
  theme_void(base_size = 8) +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1,
                                   colour = NA,
                                   size = 8),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = alpha(set1[[3]], 0.2),
                                        colour = "black")) +
  scale_y_discrete(breaks = c(1:7), expand = c(0, 0)) +
  scale_x_continuous(labels = MakeFakeLabels, expand = c(0, 0)) +
  ggtitle("Phenotype correlation")
dummy_plot

# generate a tree

p <- ggtree(phylo) +
  theme_void(base_size = 8) +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1,
                                   colour = NA,
                                   size = 8),
        plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(breaks = c(1:7)) +
  scale_x_continuous(labels = MakeFakeLabels, expand = c(0, 0.5)) +
  ggtitle("Cluster") +
  geom_tiplab()
p


layout_matrix <- matrix(
  c(rep(1, 10),
    rep(2, 28),
    rep(5, 3),
    rep(3, 28),
    rep(5, 3),
    rep(4, 28)),
  nrow = 1)
grobs <- arrangeGrob(p,
                     gp,
                     enrichment_plot,
                     dummy_plot,
                     layout_matrix = layout_matrix)


cairo_pdf("test/cluster_heatmaps.pdf",
          width = convertUnit(unit(180, "mm"), "in", valueOnly = TRUE),
          height = convertUnit(unit(90, "mm"), "in", valueOnly = TRUE))
grid.newpage()
grid.draw(grobs)
dev.off()


############################
# PLOT INDIVIDUAL CLUSTERS #
############################

# extract per-accession stats
wald_results[, accession := gsub(".*group ([[:alpha:]]+).*", "\\1", wald_test)]
acc_results <- wald_results[, .(accession, gene_id, log2FoldChange, lfcSE, symbols)]
acc_results[!symbols == "", label := paste(gene_id, symbols, sep = " | ")]
acc_results[symbols == "", label := gene_id]

# generate plot data
QuickClusterPlot <- function(cluster_number) {
  my_acc_order <- c("rufipogon",
                    "indica",
                    "japonica",
                    "barthii",
                    "glaberrima")
  file_path <- paste0("test/cluster_", cluster_number, ".pdf")
  plot_genes <- clusters[cluster == cluster_number, unique(MsuID)]
  my_fills <- RColorBrewer::brewer.pal(4, "Paired")[c(1, 2, 2, 3, 4)]
  pd2 <- acc_results[gene_id %in% plot_genes]
  pd2[, accession := factor(accession, levels = my_acc_order)]
  gp <- ggplot(pd2, aes(x = accession,
                        y = log2FoldChange,
                        ymin = log2FoldChange - lfcSE,
                        ymax = log2FoldChange + lfcSE,
                        fill = accession)) +
    theme_grey(base_size = 6) +
    xlab(NULL) + ylab("L2FC ± SE") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    scale_fill_manual(values = my_fills, guide = FALSE) +  
    facet_wrap(~ label, scales = "free_y") +
    geom_errorbar(width = 0.1) +
    geom_col(position = "dodge", alpha = 0.5)
  
  ggsave(file_path,
         gp,
         width = 10,
         height = 7.5,
         units = "in")
}

sapply(clusters[, unique(cluster)], QuickClusterPlot)



