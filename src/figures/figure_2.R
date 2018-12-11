library(data.table)
library(ggplot2)

# Figure 2: transcriptome-pca

pcx_file <- "output/050_deseq/rlog_pca/pcx.Rds"
pca_file <- "output/050_deseq/rlog_pca/pc.Rds"

long_spec_order <- c("rufipogon" = "O. rufipogon",
                     "indica" = "O. sativa indica",
                     "japonica" = "O. sativa japonica",
                     "barthii" = "O. barthii",
                     "glaberrima" = "O. glaberrima")

stage_order <- c(PBM = "IM", SM = "DM")


########
# MAIN #
########

# data
pcx <- data.table(readRDS(pcx_file))
pca <- readRDS(pca_file)

# melt 1 -> 5
pd <- melt(pcx,
     id.vars = c("accession", "stage", "continent", "domestication", "ID"),
     measure.vars = paste0("PC", 1:5),
     variable.name = "component",
     value.name = "score")

# fix order
pd[, accession := factor(plyr::revalue(accession, long_spec_order),
                          levels = long_spec_order)]
pd[, stage := factor(plyr::revalue(stage, stage_order),
                      levels = stage_order)]
pd[, rep := factor(as.numeric(gsub("[^[:digit:]]+", "", ID)))]

# calculate percent variance for facet labels
pct_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
pv_dt <- data.table(comp = paste0("PC", 1:5),
                    pv = pct_var[1:5])
pv_dt[, facet_label := paste0(comp, " (", round(pv, 1), "%)")]

pd[, facet_label := pv_dt[comp == component, facet_label], by = component]

# plot
gp <- ggplot(pd, aes(x = stage,
               y = score,
               group = ID,
               fill = accession)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(panel.background = element_rect(colour = "black"),
        strip.text.x = element_text(face= "italic")) +
  xlab(NULL) +
  ylab("Score on component") +
  facet_grid(facet_label ~ accession, scales = "free_y") +
  scale_fill_manual(values = pc5_cols[c(1, 2, 2, 3, 4)],
                    guide = FALSE) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7)

ggsave("test/Figure_2.pdf",
       gp,
       width = 178,
       height = 145,
       units = "mm")

