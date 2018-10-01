library(tidyverse)
library(readxl)
library(ggpubr)
library(gridExtra)
library(ggbeeswarm)

dat <- read_excel(path = "../data-raw/supp-table-PanicleTraitsPhenotypingAP2Mutant.xlsx")


dat <- dat %>%
  select(Id,
         `Accession Name`,
         `Pb_nb (PbN)`,
         `Sb_nb (SbN)`,
         `Sp_nb (SpN)`) %>%
  mutate_at(vars(`Pb_nb (PbN)`,
                 `Sb_nb (SbN)`),
            as.numeric)  %>%
  mutate(ID = str_split_fixed(Id, pattern = "_", n = 2)[, 1],
         locus_id = case_when(ID == "CRL5" ~ "LOC_Os07g03250",
                              ID == "SMOS1" ~ "LOC_Os05g32270",
                              TRUE ~ ""),
         is_wt = case_when(ID == "WT" ~ "WT",
                                  TRUE ~ "mutant"),
         ID = paste(ID, locus_id)) %>%
  gather(`Pb_nb (PbN)`:`Sp_nb (SpN)`,
         key = "measure",
         value = "value")


# dat <- read_excel(path = "../data-raw/Mutant_AP2_PhenotypingData.xlsx")

# gene_names <- c(LOC_Os07g03250 = "plt8",
#                 WT = "WT",
#                 LOC_Os05g32270 = "erf142",
#                 LOC_Os06g03710 = "smo2 ???",
#                 LOC_Os08g31580 = "erf48")

# dat <- dat %>% 
#   mutate_at(vars(RL:SpN), as.numeric) %>% 
#   mutate(mutant_gene = gene_names[`Target Gene`]) %>%
#   gather(RL:SpN, key = "measure", value = "value")


# Plot every mutant - Flipped -------------------------------------------------

pdf("../fig/fig-07-mutant-TEST-flipped.pdf",
    height = 4,
    width = 12)
dat %>%
  # filter(Accession != "Kitake",
  #        measure != "TbN") %>%
  # mutate(is_wt = case_when(`Target Gene` == "WT" ~ "WT",
  #                          TRUE ~ "mutant")) %>%
  ggplot(aes(x = ID,
             y = value, 
             colour = is_wt)) +
  geom_boxplot(varwidth = T,
               colour = "black",
               outlier.alpha = 0) +
  # geom_jitter(height = 0,
  #             alpha = .2) +
  geom_quasirandom(alpha = .7
                   # colour = "blue"
  ) +
  facet_grid(`Accession Name` ~ measure,
             scales = "free",
             space = "free_y") +
  coord_flip() +
  theme_bw() +
  # theme(axis.text.x = element_text(hjust = 1,
  #                                  vjust = .5)) +
  scale_color_viridis_d(begin = .2, end = .8) +
  labs(x = "Value",
       y = "Accession_mutant",
       caption = "")
# ylim(0, NA)

# ggplot(aes(x = ID,
#            y = value, 
#            colour = is_wt)) +
#   geom_boxplot(varwidth = T,
#                colour = "black",
#                outlier.alpha = 0) +
#   # geom_jitter(height = 0,
#   #             alpha = .2) +
#   geom_quasirandom(alpha = .7
#                    # colour = "blue"
#                    ) +
#   facet_grid(Accession ~ measure,
#              scales = "free",
#              space = "free_y") +
#   coord_flip() +
#   theme_bw() +
#   # theme(axis.text.x = element_text(hjust = 1,
#   #                                  vjust = .5)) +
#   scale_color_viridis_d(begin = .2, end = .8) +
#   labs(x = "Value",
#        y = "Accession_mutant",
#        caption = "")
#   # ylim(0, NA)
dev.off()


# Plot every mutant - standard --------------------------------------------

pdf("../fig/fig-07-mutant-TEST.pdf",
    height = 12,
    width = 5)
dat %>%
  # filter(Accession != "Kitake",
  #        measure != "TbN") %>%
  # mutate(is_wt = case_when(`Target Gene` == "WT" ~ "WT",
  #                          TRUE ~ "mutant")) %>%
  ggplot(aes(x = ID,
             y = value, 
             colour = is_wt)) +
  geom_boxplot(varwidth = T,
               colour = "black",
               outlier.alpha = 0) +
  # geom_jitter(height = 0,
  #             alpha = .2) +
  geom_quasirandom(alpha = .7
                   # colour = "blue"
  ) +
  facet_grid(measure ~ `Accession Name`,
             scales = "free",
             space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270,
                                   hjust = 0,
                                   vjust = .5)) +
  scale_color_viridis_d(begin = .2, end = .8) +
  labs(y = "Value",
       x = "Accession_mutant",
       caption = "")
# ylim(0, NA)


# dat %>%
#   filter(Accession != "Kitake",
#          measure != "TbN") %>%
#   mutate(is_wt = case_when(`Target Gene` == "WT" ~ "WT",
#                            TRUE ~ "mutant")) %>%
#   ggplot(aes(x = ID,
#              y = value, 
#              colour = is_wt)) +
#   geom_boxplot(varwidth = T,
#                colour = "black",
#                outlier.alpha = 0) +
#   # geom_jitter(height = 0,
#   #             alpha = .2) +
#   geom_quasirandom(alpha = .7
#                    # colour = "blue"
#   ) +
#   facet_grid(measure ~ Accession,
#              scales = "free",
#              space = "free_x") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 270,
#                                    hjust = 0,
#                                    vjust = .5)) +
#   scale_color_viridis_d(begin = .2, end = .8) +
#   labs(y = "Value",
#        x = "Accession_mutant",
#        caption = "")
dev.off()

# Plot every mutant - selected --------------------------------------------

pdf("../fig/fig-07-mutant-TEST-reduced.pdf",
    height = 7,
    width = 6)
dat %>%
  filter(Accession != "Kitake",
         measure %in% c("PbN", "SbN", "SpN")) %>%
  mutate(is_wt = case_when(`Target Gene` == "WT" ~ "WT",
                           TRUE ~ "mutant"),
         mutant_gene = case_when(mutant_gene != "WT" ~ paste(mutant_gene, 
                                                             `Target Gene`,
                                                             sep = "\n"),
                                 TRUE ~ "WT")) %>%
  ggplot(aes(x = mutant_gene,
             y = value, 
             colour = is_wt)) +
  geom_boxplot(varwidth = T,
               colour = "black",
               outlier.alpha = 0) +
  # geom_jitter(height = 0,
  #             alpha = .2) +
  geom_quasirandom(alpha = .7
                   # colour = "blue"
  ) +
  facet_grid(measure ~ Accession,
             scales = "free",
             space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270,
                                   hjust = 0,
                                   vjust = .5)) +
  scale_color_viridis_d(begin = .2, end = .8) +
  labs(x = "Genotype",
       y = "Counts [n]",
       caption = str_wrap("Mutants of four AP2-EREBP genes have panicle
                          phenotypes. When compared to the wild type,
                          erf48 produces slightly more primary
                          branches.
                          Instead, plt8 produces less primary branches and 
                          slightly less spikelets. The mutants of the two 
                          homologs smos1 (erf142) and smos2 both produce less primary
                          and secondary branches and less spikelets (and
                          overall smaller panicles - not shown).",
                          width = 80),
       colour = "")
# ylim(0, NA)
dev.off()



# Plot expression of selected genes ---------------------------------------

source("helper-functions.R")
dds <- readRDS("../data-raw/dds.Rds")

dat %>%
  filter(dat$`Target Gene` != "WT") %>%
  pull(`Target Gene`) %>%
  unique() %>%
  get_expression(dds = dds) %>%
  plot_norm_expr()
  
# dat %>%
#   filter(dat$`Target Gene` != "WT") %>%
#   select(`Target Gene`, mutant_gene) %>%
#   distinct()

# How to add pvalue to a plot? (ggpubr)--------------------------------------

tst <- dat %>%
  split(.$Accession)

tst <- dat %>%
  filter(Accession != "Kitake") %>%
  arrange(desc(mutant_gene)) %>%
  mutate(mutant_gene = as_factor(mutant_gene)) %>%
  split(paste(.$Accession)) %>%
  map(~split(., .$measure))

p <- list()
p$Nipponbare <- map(tst$Niponbarre, ~ggboxplot(., x = "mutant_gene", y = "value",
                                               ylab = unique(.$measure),
                                               title = unique(.$measure),
                                               fill = "cornsilk") +
                      stat_compare_means(method = "t.test",
                                         comparisons = list(c(1,2), c(2,3), c(1,3))))

p$Kinmaze <- map(tst$Kinmaze, ~ggboxplot(., x = "mutant_gene", y = "value",
                                         ylab = unique(.$measure),
                                         title = unique(.$measure),
                                         fill = "cornsilk") +
                   stat_compare_means(method = "t.test",
                                      comparisons = list(c(1,2))))

p$Illmibyeo <- map(tst$Illmibyeo, ~ggboxplot(., x = "mutant_gene", y = "value",
                                             ylab = unique(.$measure),
                                             title = unique(.$measure),
                                             fill = "cornsilk") +
                     stat_compare_means(method = "t.test",
                                        comparisons = list(c(1,2))))

arrange_plots <- function(name) {
  header <- text_grob(name,
                      face = "bold",
                      size = 32)
  p1 <- ggarrange(plotlist = p[[name]],
                  ncol = length(p[[name]])/3,
                  nrow = length(p[[name]])/3)
  ggarrange(header, p1, nrow = 2, 
            heights = c(1, 9))
}

pdf(file = "../fig/fig-07-mutants_with_stats.pdf",
    width = 12,
    height = 14)
map(names(p), arrange_plots)
dev.off()
