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

set.seed(1)

phenotypes_file <- snakemake@input[["phenotypes"]]
pheno_names_file <- snakemake@input[["pheno_names"]]

# plots
fig1_file <- snakemake@output[["fig1"]]

# dev
# phenotypes_file <- "data/phenotyping/panicle_mutants.csv"
# pheno_names_file <- "data/phenotyping/phenotype_name_key.csv"

########
# MAIN #
########

# data
phenotypes <- fread(phenotypes_file)
pheno_names <- fread(pheno_names_file)

# tidy
old_names <- names(phenotypes)[names(phenotypes) %in% pheno_names[, mtp_name]]
new_names <- pheno_names[old_names, on = .(mtp_name), text_name]
setnames(phenotypes,
         old_names,
         new_names)

setnames(phenotypes,
         c("Id",
           "Plant_nb (1 to 3)",
           "Panicle_nb (1 to 3)",
           "Type (WT/mutant)",
           "Species",
           "Accession Name"),
         c("id",
           "plant_number",
           "panicle_number",
           "type",
           "species",
           "accession"))

pheno <- melt(phenotypes,
              id.vars = c("id",
                          "plant_number",
                          "panicle_number",
                          "accession"),
              measure.vars = c("Primary branch number",
                               "Secondary branch number",
                               "Spikelet number"))

# fix genotype
pheno[, gt := gsub("_.*", "", id)]
pheno[gt != "WT", gt := paste0('italic("', tolower(gt), '")')]
pheno[gt == "WT", gt := '"WT"']

# arrange colour
pheno[, pt_col := gt == '"WT"']

#draw the plot
gp <- ggplot(pheno, aes(x = gt, y = value, colour = pt_col)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(panel.background = element_rect(colour = "black"),
        strip.placement = "outside") +
  facet_grid(variable ~ accession, scales = "free", switch = "y") +
  xlab(NULL) + ylab(NULL) +
  scale_x_discrete(labels = function(l) parse(text = l)) +
  scale_color_brewer(palette = "Set1",
                     guide = FALSE) +
  geom_boxplot(colour = alpha("black", 0.5),
               fill = NA,
               outlier.colour = NA,
               width = 0.55) +
  geom_point(position = position_jitter(width = 0.3, height = 0),
             size = 1.5,
             alpha = 0.8,
             shape = 16) 

ggsave(fig1_file,
       gp,
       width = 87,
       height = 87*3/2,
       units = "mm")

# Log
sessionInfo()
