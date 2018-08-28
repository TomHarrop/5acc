library(cowplot) # add to MS repo!!!
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(png)

names_file <- "data/phenotyping/phenotype_name_key.csv"
pca_file <- "output/080_phenotype/cali_pca.Rds"
pheno_file <- "output/080_phenotype/cali.csv"

########
# MAIN #
########

# read data
pca <- readRDS(pca_file)
pheno <- fread(pheno_file)
pheno_names <- fread(names_file)
b_png <- readPNG("test/test_resized.png")

# make plotting data
spec_order <- c("rufipogon" = "O. rufipogon",
                "indica" = "O. sativa",
                "japonica" = "O. sativa",
                "sativa" = "O. sativa",
                "barthii" = "O. barthii",
                "glaberrima" = "O. glaberrima")


pca_pheno <- cbind(pheno, pca$x)
pca_pd <- melt(pca_pheno,
               id.vars = c("Species",
                           "Sowing_nb",
                           "Repet_nb",
                           "Plant_nb",
                           "Panicle_nb"),
               measure.vars = paste0("PC", 1:4),
               variable.name = "component")

pca_pd[, Species := factor(plyr::revalue(Species, spec_order),
                           levels = unique(spec_order))]

# loadings
loadings_wide <- data.table(pca$rotation, keep.rownames = TRUE)
setnames(loadings_wide, "rn", "phenotype")

# cluster loadings for plot order
hc <- hclust(dist(pca$rotation, method = "minkowski"),
             method = "ward.D2")
pheno_order <- rev(hc$labels[hc$order])
setkey(pheno_names, full_name)
names_order <- pheno_names[pheno_order, short_name]

# generate loadings plot data
loadings_pd <- melt(loadings_wide,
                    id.vars = "phenotype",
                    measure.vars = paste0("PC", 1:4),
                    variable.name = "component")
loadings_pd[, phenotype := factor(phenotype, levels = pheno_order)]
loadings_pd[, phenotype := plyr::mapvalues(phenotype,
                                           pheno_order,
                                           names_order)]

# plot pca
pd <- RColorBrewer::brewer.pal(4, "Paired")

pcp <- ggplot(pca_pd, aes(x = Species, y = value, colour = Species)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = "italic"),
        panel.border = element_rect(fill = NA, colour = "black")) +
  facet_wrap(~ component, nrow = 1) +
  xlab(NULL) +
  ylab("Value") +
  #ggtitle("PCA on Cali phenotyping", label = "(A)") +
  scale_color_manual(values = pd[c(1, 2, 3, 4)],
                     guide = FALSE) +
  geom_point(position = position_jitter(width = 0.4),
             size = 1,
             alpha = 0.8,
             shape = 16) +
  geom_boxplot(colour = alpha("black", 0.5),
               fill = NA,
               outlier.colour = NA)
pcp

# plot loadings
lp <- ggplot(loadings_pd,
             aes(x = phenotype,
                 yend = value,
                 colour = phenotype,
                 xend = phenotype)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(fill = NA, colour = "black")) +
  xlab(NULL) + ylab("Loading") +
  #ggtitle(label = "(B)") +
  scale_color_brewer(palette = "Set1", guide = FALSE) +
  facet_wrap(~ component, nrow = 1) +
  geom_segment(y = 0,
               size = 1,
               arrow = arrow(angle = 15, length = unit(0.1, "in")))
lp


# cowplot
vp <- viewport(width = unit(1, "npc"), height = unit(0.95, "npc"))
b_grob <- rasterGrob(b_png, interpolate = TRUE) 
barthi_grob <- ggplot(data.frame(x = 0.5, y = 0, label = 'italic("O. barthii")'),
       aes(x = x, y = y, label = label)) +
  theme_void(base_size = 8, base_family = "Helvetica") +
  theme(plot.margin = unit(c(0, 0, 2, 0), "mm")) +
  xlim(c(0, 1)) + ylim(c(0, 1)) +
  annotation_custom(b_grob) +
  geom_text(parse = TRUE)

lhs <- plot_grid(barthi_grob, barthi_grob, barthi_grob, barthi_grob,
          labels = c("A", "B", "C", "D"),
          ncol = 1,
          label_size = 10,
          label_fontfamily = "Helvetica")

pca_plots <- plot_grid(pcp, lp,
                     ncol = 1,
                     align = "hv",
                     axis = "tblr",
                     labels = c("E", "F"),
                     label_size = 10,
                     label_fontfamily = "Helvetica")

cowplot <- plot_grid(lhs, pca_plots, ncol = 2, rel_widths = c(1, 2))

ggsave("test/pc1-4_cowplot.pdf", cowplot, width = 178, height = 225, units = "mm")


# draw plots
layout_matrix <- matrix(
  1:2,
  nrow = 2)

grobs <- arrangeGrob(pcp,
                     lp,
                     layout_matrix = layout_matrix)

# save output
cairo_pdf("test/pc1-4.pdf",
          width = convertUnit(unit(180, "mm"), "in", valueOnly = TRUE),
          height = convertUnit(unit(180, "mm"), "in", valueOnly = TRUE))
grid.newpage()
grid.draw(grobs)
dev.off()
