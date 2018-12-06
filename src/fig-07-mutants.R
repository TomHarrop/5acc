library(tidyverse)
library(readxl)
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
         ID = case_when(ID == "CRL5" ~ tolower(ID),
                        ID == "SMOS1" ~ tolower(ID),
                        TRUE ~ ID),
         is_wt = case_when(ID == "WT" ~ "WT",
                                  TRUE ~ "mutant")
         # ID = paste(tolower(ID), locus_id, sep = "\n"),
         ) %>%
  gather(`Pb_nb (PbN)`:`Sp_nb (SpN)`,
         key = "measure",
         value = "value") %>%
  mutate(measure = case_when(measure == "Pb_nb (PbN)" ~ "Primary Branches",
                             measure == "Sb_nb (SbN)" ~ "Secondary Branches",
                             measure == "Sp_nb (SpN)" ~ "Spikelets",
                             TRUE ~ measure),
         `Accession Name` = case_when(`Accession Name` == "Niponbarre" ~ "Nipponbare",
                                      TRUE ~ `Accession Name`)) %>%
  ## Mutant first in X axis
  arrange(desc(ID)) %>%
  mutate(ID = as_factor(ID))


# Plot every mutant - standard --------------------------------------------

p <-
  dat %>%
  ggplot(aes(x = ID,
             y = value, 
             colour = is_wt)) +
  geom_boxplot(varwidth = T,
               colour = "black",
               outlier.alpha = 0) +
  geom_quasirandom(size = 2,
                   alpha = .8) +
  facet_grid(measure ~ `Accession Name`,
             scales = "free",
             space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270,
                                   hjust = 0,
                                   vjust = .5,
                                   face = "italic"), 
        legend.position = "top",
        legend.background = element_rect(size=0.2,
                                         linetype="solid",
                                         colour = "grey80")) +
  scale_color_viridis_d(begin = .2,
                        end = .8,
                        guide = FALSE) +
  # guides(colour = guide_legend(title = NULL,
  #                             label.position = "right",
  #                             nrow=1,
  #                             override.aes = list(alpha = 1))) +
  labs(# title = "Phenotype of AP2 Mutants",
       # caption = str_wrap("Mutants of two AP2/EREBP-like genes,
       #                    CRL5 and SMOS1, have defects in panicle architecture.
       #                    The crl5 mutant produces fewer primary branches.
       #                    The smos1 mutant produces fewer primary branches,
       #                    secondary branches and spikelets.
       #                    These two mutants are in different O. sativa genetic background:
       #                    crl5 is in Kinmaze background and smos1 is in Nipponbare background.
       #                    They are both compared with their own background WT variety.",
       #                    width = 50),
       y = "Count (n)",
       x = "")

pdf("../fig/fig- panicle-mutants.pdf",
    height = 7,
    width = 3.4)
p %>% print()
dev.off()


# Plot every mutant - Flipped -------------------------------------------------

# p_flipped <- 
#   dat %>%
#   ggplot(aes(x = ID,
#              y = value, 
#              colour = is_wt)) +
#   geom_boxplot(varwidth = T,
#                colour = "black",
#                outlier.alpha = 0) +
#   geom_quasirandom(alpha = .7) +
#   facet_grid(`Accession Name` ~ measure,
#              scales = "free",
#              space = "free_y") +
#   coord_flip() +
#   theme_bw() +
#   scale_color_viridis_d(begin = .2, end = .8) +
#   labs(x = "Value",
#        y = "Accession_mutant",
#        caption = "")
# 
# 
# pdf("../fig/fig-07-mutant-TEST-flipped.pdf",
#     height = 4,
#     width = 12)
# p_flipped %>% print()
# dev.off()

