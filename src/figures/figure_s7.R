#!/usr/bin/env Rscript

library(DESeq2)
library(data.table)
library(ggplot2)
library(viridisLite)

# Figure S7: Heatmap of pairwise distances between libraries.

dds_file <- "output/050_deseq/filtered_dds.Rds"

spec_order <- c("or" = "O. rufipogon",
                "osi" = "O. sativa indica",
                "osj" = "O. sativa japonica",
                "ob" = "O. barthii",
                "og" = "O. glaberrima")

stage_order <- c(PBM = "IM", SM = "DM")

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

# split labels
SampleToLabel <- function(x){
  my_split <- unlist(strsplit(x, "_"))
  my_label <- paste0('italic("',
                     spec_order[my_split[1]],
                     '")~"',
                     stage_order[my_split[2]],
                     '"~"',
                     my_split[3],
                     '"')
  return(my_label)
}
dist_dt[, l1_label := SampleToLabel(as.character(l1)), by = l1]
dist_dt[, l2_label := SampleToLabel(as.character(l2)), by = l2]

# sample order
hc <- hclust(sample_dist, method = "ward.D2")
sample_order <- hc$labels[hc$order]
dist_dt[, l1 := factor(l1, levels = sample_order)]
dist_dt[, l2 := factor(l2, levels = sample_order)]
setorder(dist_dt, l1, l2)
dist_dt[, l1_label := factor(l1_label, levels = unique(l1_label))]
dist_dt[, l2_label := factor(l2_label, levels = unique(l2_label))]

# plot
gp <- ggplot(dist_dt, aes(x = l1_label, y = l2_label, fill = Distance)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.length = unit(0, "mm")) +
  coord_fixed() +
  scale_x_discrete(expand = c(0,0),
                   labels = function(l) parse(text = l)) +
  scale_y_discrete(expand = c(0,0),
                   labels = function(l) parse(text = l)) +
  xlab(NULL) + ylab(NULL) +
  scale_fill_viridis_c() +
  geom_raster()

ggsave("test/Figure_S7.pdf",
       device = cairo_pdf,
       gp,
       width = 178,
       height = 178,
       units = "mm")

