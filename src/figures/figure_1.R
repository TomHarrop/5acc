#!/usr/bin/env Rscript

# Figure 1: panicle photos with Cali PCA

library(cowplot) # add to MS repo!!!
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(png) # add to MS repo!!!

#############
# FUNCTIONS #
#############

PlotPng <- function(my_png, my_label){
  
  my_grob <- rasterGrob(my_png, interpolate = TRUE) 
  my_expr <- glue::glue('italic("{my_label}")')
  
  my_plot <- ggplot(data.frame(x = 1, y = 1, label = my_expr),
                    aes(x = x, y = y, label = label)) +
    theme_void(base_size = 8, base_family = "Helvetica") +
    theme(plot.margin = unit(c(0, 1, 0, 1), "mm")) + #trbl
    xlim(c(0, 1)) + ylim(c(0, 1)) +
    annotation_custom(my_grob) +
    geom_text(parse = TRUE,
              size = 2.5,
              hjust = "inward",
              nudge_x = 0)
  
  return(my_plot)
}


###########
# GLOBALS #
###########

names_file <- "data/phenotyping/phenotype_name_key.csv"
pca_file <- "output/080_phenotype/cali_pca.Rds"
pheno_file <- "output/080_phenotype/cali.csv"

# png files
ob_pan_file <- "test/panicles/Ob_B88.png"
os_pan_file <- "test/panicles/Os_IR64.png"
og_pan_file <- "test/panicles/Og_Tog5681.png"
or_pan_file <- "test/panicles/Or_W1654.png"

########
# MAIN #
########

# read data
pca <- readRDS(pca_file)
pheno <- fread(pheno_file)
pheno_names <- fread(names_file)

# pngs
ob_png <- readPNG(ob_pan_file)
os_png <- readPNG(os_pan_file)
or_png <- readPNG(or_pan_file)
og_png <- readPNG(og_pan_file)

# calculate percent variance for facet labels
pct_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
pv_dt <- data.table(component = paste0("PC", 1:9),
           pv = pct_var[1:9])
pv_dt[, facet_label := paste0(component, " (", round(pv, 1), "%)")]

# make plotting data
spec_order <- c("rufipogon" = "O. rufipogon",
                "indica" = "O. sativa",
                "japonica" = "O. sativa",
                "sativa" = "O. sativa",
                "barthii" = "O. barthii",
                "glaberrima" = "O. glaberrima")

pca_pheno <- cbind(pheno, pca$x)
pca_pd_long <- melt(pca_pheno,
               id.vars = c("Species",
                           "Sowing_nb",
                           "Repet_nb",
                           "Plant_nb",
                           "Panicle_nb"),
               measure.vars = paste0("PC", 1:9),
               variable.name = "component")

pca_pd <- merge(pca_pd_long, pv_dt)
pca_pd[, Species := factor(plyr::revalue(Species, spec_order),
                           levels = unique(spec_order))]

# loadings
loadings_wide <- data.table(pca$rotation, keep.rownames = TRUE)
setnames(loadings_wide, "rn", "phenotype")

# order loadings by contribution to PC1
order_dt <- data.table(pca$rotation, keep.rownames = TRUE)
setorder(order_dt, -PC1)
pheno_order <- order_dt$rn
setkey(pheno_names, full_name)
names_order <- pheno_names[pheno_order, short_name]

# generate loadings plot data
loadings_pd <- melt(loadings_wide,
                    id.vars = "phenotype",
                    measure.vars = paste0("PC", 1:9),
                    variable.name = "component")
loadings_pd[, phenotype := factor(phenotype, levels = pheno_order)]
loadings_pd[, phenotype := plyr::mapvalues(phenotype,
                                           pheno_order,
                                           names_order)]

# plot pca
pd <- RColorBrewer::brewer.pal(4, "Paired")

pcp1 <- ggplot(pca_pd[component == "PC1"],
              aes(x = Species, y = value, colour = Species)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        panel.border = element_rect(fill = NA, colour = "black")) +
  facet_wrap(~ facet_label, nrow = 1) +
  xlab(NULL) +
  ylab("Value") +
  scale_color_manual(values = pd[c(1, 2, 3, 4)],
                     guide = FALSE) +
  geom_point(position = position_jitter(width = 0.4),
             size = 1,
             alpha = 0.8,
             shape = 16) +
  geom_boxplot(colour = alpha("black", 0.5),
               fill = NA,
               outlier.colour = NA)
pcp_all <- ggplot(pca_pd[pv > 10],
               aes(x = Species, y = value, colour = Species)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        panel.border = element_rect(fill = NA, colour = "black")) +
  facet_wrap(~ facet_label, nrow = 1) +
  xlab(NULL) +
  ylab("Value") +
  scale_color_manual(values = pd[c(1, 2, 3, 4)],
                     guide = FALSE) +
  geom_point(position = position_jitter(width = 0.4),
             size = 1,
             alpha = 0.8,
             shape = 16) +
  geom_boxplot(colour = alpha("black", 0.5),
               fill = NA,
               outlier.colour = NA)
pcp1
pcp_all

# plot loadings
lp1 <- ggplot(loadings_pd[component == "PC1"],
             aes(x = phenotype,
                 yend = value,
                 colour = phenotype,
                 xend = phenotype)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(fill = NA, colour = "black")) +
  xlab(NULL) + ylab("Loading") +
  scale_color_brewer(palette = "Set1", guide = FALSE) +
  facet_wrap(~ component, nrow = 1) +
  geom_segment(y = 0,
               size = 1,
               arrow = arrow(angle = 15, length = unit(0.1, "in")))
lp_all <- ggplot(loadings_pd[component %in% pca_pd[pv > 10, unique(component)]],
              aes(x = phenotype,
                  yend = value,
                  colour = phenotype,
                  xend = phenotype)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(fill = NA, colour = "black")) +
  xlab(NULL) + ylab("Loading") +
  scale_color_brewer(palette = "Set1", guide = FALSE) +
  facet_wrap(~ component, nrow = 1) +
  geom_segment(y = 0,
               size = 1,
               arrow = arrow(angle = 15, length = unit(0.1, "in")))

lp1
lp_all

# plot pngs
ob_plot <- PlotPng(ob_png, "O. barthii")
os_plot <- PlotPng(os_png, "O. sativa")
or_plot <- PlotPng(or_png, "O. rufipogon")
og_plot <- PlotPng(og_png, "O. glaberrima")


# layout figure 1
top <- plot_grid(or_plot, os_plot, ob_plot, og_plot,
                 labels = c("A", "B", "C", "D"),
                 align = "hv",
                 axis = "tblr",
                 ncol = 4,
                 label_size = 10,
                 label_fontfamily = "Helvetica")

bottom <- plot_grid(pcp1, lp1,
                    ncol = 2,
                    align = "hv",
                    axis = "tblr",
                    labels = c("E", "F"),
                    label_size = 10,
                    label_fontfamily = "Helvetica",
                    rel_widths = c(3, 1))

cowplot <- plot_grid(top, bottom, align = "v",
                     axis = "l",
                     
          ncol = 1, 
          rel_heights = c(1, 1))

ggsave("test/Figure_1.pdf",
       cowplot,
       width = 178,
       height = 225/2,
       units = "mm")

# layout SF1
cowplot2 <- plot_grid(pcp_all,
                      lp_all,
                      nrow = 2,
                      align = "hv",
                      axis = "tlbr")
ggsave("test/SF1.pdf",
       cowplot2,
       width = 178,
       height = 225/2,
       units = "mm")


