
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