#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)


phenotypes_file <- "data/phenotyping/Phenotype_PanicleSequenced_corrected2.csv"
pheno_names_file <- "data/phenotyping/phenotype_name_key.csv"

# data
pheno_wide <- fread(phenotypes_file)
pheno_names <- fread(pheno_names_file)

# tidy
setnames(pheno_wide,
         pheno_names[mtp_name != "", mtp_name],
         pheno_names[mtp_name != "", full_name])

spec_order <- c("Oryza rufipogon" = "O. rufipogon",
                "Oryza sativa indica" = "O. sativa indica",
                "Oryza sativa japonica temp" = "O. sativa japonica",
                "Ozyza Barthii" = "O. barthii",
                "Oryza glaberrima" = "O. glaberrima")


pheno <- melt(pheno_wide,
              id.vars = "Species",
              measure.vars = c("primary_branch_number",
                               "secondary_branch_number",
                               "spikelet_number",
                               "primary_branch_length",
                               "secondary_branch_length",
                               "rachis_length"))

pheno[, species_label := factor(plyr::revalue(Species, spec_order),
                                levels = spec_order)]

pheno_pd <- merge(pheno,
                  pheno_names,
                  by.x = "variable",
                  by.y = "full_name",
                  all.x = TRUE,
                  all.y = FALSE)
panel_order <- c("Rachis length",
                 "Primary branch number",
                 "Primary branch length", 
                 "Secondary branch number",
                 "Secondary branch length",
                 "Spikelet number")
panel_names <- c("Rachis length (cm)",
                 "Primary branch number",
                 "Primary branch length (cm)", 
                 "Secondary branch number",
                 "Secondary branch length (cm)",
                 "Spikelet number")
pheno_pd[, text_name := factor(plyr::mapvalues(text_name,
                                               panel_order,
                                               panel_names),
                               levels = panel_names)]

# draw plot
paired <- RColorBrewer::brewer.pal(4, "Paired")
gp <- ggplot(pheno_pd, aes(x = species_label, y = value, colour = species_label)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(panel.background = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   face = "italic"),
        strip.placement = "outside") +
  scale_color_manual(values = paired[c(1, 2, 2, 3, 4)],
                     guide = FALSE) +
  xlab(NULL) + ylab(NULL ) +
  facet_wrap(~text_name, scales = "free_y", strip.position = "left", nrow = 3) +
  geom_point(position = position_jitter(width = 0.3),
             size = 1,
             alpha = 0.8,
             shape = 16)

ggsave("test/Figure_S5.pdf",
       device = cairo_pdf,
       gp,
       width = 87,
       height = 130,
       units = "mm")

