library(tidyverse)
library(ggfortify)

load("../data/phenotypes.Rdata")


# How many plants and how many panicles -----------------------------------

pheno_cali %>%
  group_by(Bar_Code) %>%
  summarise(plants = max(Plant_nb),
            panicles = max(Panicle_nb))

# PCA ---------------------------------------------------------------------

pc <- pheno_cali %>%
  select_at(vars(RL:spn)) %>%
  select(-PanL) %>%
  as.data.frame() %>%
  prcomp(.,
         scale. = T,
         center = T)

p <- autoplot(pc, data = pheno_cali, 
              loadings = TRUE,
              loadings.label = TRUE) 

# Plot boxplot by accession -----------------------------------------------

p2 <- 
  p$data %>%
  mutate(Name = case_when(Name == "Nipponbare" ~ "Niponbarre",
                          Name == "W1654 (2) / B" ~ "W1654",
                          TRUE ~ Name),
         in_rnaseq = case_when(Name %in%
                                 pheno_mnp$species ~ "yes",
                               TRUE ~ "no"),
         Origin = case_when(Origin == "Ob" ~ "O. barthii",
                            Origin == "Og" ~ "O. glaberrima",
                            Origin == "Or" ~ "O. rufipogon",
                            Origin == "Os" ~ "O. sativa")) %>%
  ggplot(aes(x = reorder(Name, PC1), 
             y = PC1,
             fill = in_rnaseq)) +
  geom_boxplot() +
  facet_grid(. ~ Origin,
             scales = "free_x",
             space = "free_x") +
  scale_fill_manual(values = c("white", "red")) +
  # coord_flip() +
  labs(title = paste0("Loadings of PC1 for all varieties of rice ",
                        "that we phenotyped"), 
       x = "Rice Accession",
       y = "Loading on PC1",
       fill = "In RNA-seq", 
       caption = str_wrap("Principal component analysis (PCA) of panicle phenotyping.
                          PC1 splits wild and domesticated varieties.
                          The varieties included in the RNAseq experiment are displayed 
                          in red. Those varieties tend to average the characterstic
                          of their species of origin (all but Nipponbare, which was
                          included because as reference). 
                          For each variety, we have measured from 3 to 9 panicles
                          (3 panicles from each of 1 to 3 plants).")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = .5,
                                   angle = 270),
        strip.text = element_text(face = "italic")) 

# a4 = 8.27 × 11.69 inches
pdf(file = "../fig/suppl-fig-phenotype-all-varieties.pdf",
    width = 11,
    height = 5.7,
    paper = "a4r")
p2 %>% print()
dev.off()


# OLD ---------------------------------------------------------------------

# PCA with labels ----------------------------------------------------------

# pc <- pheno_cali %>%
#   select_at(vars(Name, RL:spn)) %>%
#   # select(-PanL) %>%
#   as.data.frame() %>%
#   group_by(Name) %>%
#   mutate(cnt = 1:n()) %>%
#   ungroup() %>%
#   mutate(cnt = paste(Name, cnt, sep = "-")) %>%
#   as.data.frame() %>%
#   column_to_rownames("cnt") %>%
#   select(-Name) %>%
#   prcomp(.,
#          scale. = T,
#          center = T)
# 
# p <- autoplot(pc, data = pheno_cali, 
#               loadings = TRUE,
#               loadings.label = TRUE) 
# 
# pos <- position_jitter(width = .4,
#                        height = 0, 
#                        seed = 1)
# 
# 
# # Plot PCA with labels ----------------------------------------------------
# 
# 
# 
# p$data %>%
#   mutate(sequenced = case_when(Name %in% pheno_mnp$species ~ "in rnaseq",
#                                TRUE ~ "no")) %>%
#   ggplot(aes(x = Origin, 
#              y = PC1)) +
#   geom_boxplot(outlier.alpha = 0, colour = "black") +
#   geom_jitter(data = . %>%
#                 filter(sequenced == "no"),
#               height = 0,
#               alpha = .5,
#               colour = "darkgrey") +
#   geom_point(data = . %>%
#                filter(sequenced == "in rnaseq"),
#              # alpha = .5,
#              colour = "red",
#              position = pos) + 
#   # scale_color_manual(values = c("red", "grey")) +
#   ggrepel::geom_label_repel(data = . %>%
#                               filter(sequenced == "in rnaseq"),
#                             aes(label = Name),
#                             position = pos) +
#   theme_bw()
# 
# 
# 
