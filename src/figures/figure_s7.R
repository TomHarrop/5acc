#!/usr/bin/env Rscript

library(DESeq2)
library(data.table)
library(ggplot2)
library(viridisLite)

# Figure S7: Heatmap of pairwise distances between libraries.

dds_file <- "output/050_deseq/filtered_dds.Rds"

# transform read counts
dds <- readRDS(dds_file)
vst <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_assay <- assay(vst)

# pairwise sample distances
sample_dist <- dist(t(vst_assay), method = "minkowski")
dist_matrix <- as.matrix(sample_dist)
sample_dists_wide <- data.table(dist_matrix,   keep.rownames = TRUE)

# long dt
dist_dt <- melt(sample_dists_wide,
     id.vars = "rn",
     variable.name = "l2",
     value.name = "Distance")
setnames(dist_dt, "rn", "l1")
dist_dt <- unique(dist_dt, by = c("l1", "l2"))

# sample order
hc <- hclust(sample_dist, method = "ward.D2")
sample_order <- hc$labels[hc$order]
dist_dt[, l1 := factor(l1, levels = sample_order)]
dist_dt[, l2 := factor(l2, levels = sample_order)]
setorder(dist_dt, l1, l2)

# plot
gp <- ggplot(dist_dt, aes(x = l1, y = l2, fill = Distance)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  scale_fill_viridis_c() +
  geom_raster()

ggsave("test/Figure_S7.pdf",
       device = cairo_pdf,
       gp,
       width = 178,
       height = 100,
       units = "mm")

