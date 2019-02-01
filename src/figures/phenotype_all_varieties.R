#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)

###########
# GLOBALS #
###########

pheno_cali_file <- snakemake@input[["pheno_cali"]]
cali_pca_file <- snakemake@input[["cali_pca"]]

# plots
sf1_file <- snakemake@output[["sf1"]]

spec_order <- c("rufipogon" = "O. rufipogon",
                "indica" = "O. sativa",
                "japonica" = "O. sativa",
                "sativa" = "O. sativa",
                "barthii" = "O. barthii",
                "glaberrima" = "O. glaberrima")

# dev
# pheno_cali <- fread("output/080_phenotype/cali.csv")
# cali_pca <- readRDS("output/080_phenotype/cali_pca.Rds")

########
# MAIN #
########

pheno_cali <- fread(pheno_cali_file)
cali_pca <- readRDS(cali_pca_file)

# combine PCA with names
pca_pheno <- cbind(pheno_cali, cali_pca$x)

# order by PC1
mean_pca <- pca_pheno[, median(PC1), by = Name]
setorder(mean_pca, V1)
name_order <- mean_pca[, as.character(unique(Name))]
pca_pheno[, Name := factor(Name, levels = name_order)]

# order species
pca_pheno[, Species := factor(plyr::revalue(Species, spec_order),
                           levels = unique(spec_order))]

# paint some red
in_rnaseq <- c("Nipponbare", "IR64", "W1654", "Tog5681", "B88")
pca_pheno[, in_rnaseq := Name %in% in_rnaseq]

# plot
gp <- ggplot(pca_pheno, aes(x = Name, y = PC1, fill = in_rnaseq)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(strip.text = element_text(face = "italic"),
        panel.background = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 6)) +
  xlab(NULL) + ylab("Value on PC1") +
  scale_fill_manual(values = c(NA, "red"),
                    guide = FALSE) +
  facet_grid(. ~ Species, scales = "free_x", space = "free_x") +
  geom_boxplot(outlier.size = 1)

ggsave(sf1_file,
       gp,
       device = cairo_pdf,
       width = 178,
       height = 225/2,
       units = "mm")

# Log
sessionInfo()


