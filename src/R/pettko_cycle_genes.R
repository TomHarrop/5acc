#!/usr/bin/env Rscript

library(data.table)
source("src/R/column_plots.R")

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}
GenerateMessage("Plotting pettko cycle genes")

# load file
pettko.file <- "data/goi/pettko_cycle.genes.Rds"
if (!file.exists(pettko.file)) {
  stop("Couldn't find pettko_cycle.genes.Rds")
}
pettko.genes <- readRDS(pettko.file)

# get tpm values
tpm.file <- "output/tpm/tpm_with_calls.Rds"
if (!file.exists(tpm.file)) {
  stop("Couldn't find tpm.Rds")
}
tpm.long <- readRDS(tpm.file)

ord <- c("O. rufipogon", "O. sativa indica", "O. sativa japonica",
         "O. barthii", "O. glaberrima")
tpm.long[, accession := substr(sample, 1, 1)]
tpm.long[ , accession := factor(plyr::mapvalues(
  accession,
  from = c("R", "I", "J", "B", "G"),
  to = ord),  levels = ord)
  ]

###############
# COLUMN PLOT #
###############

# generate data for plotting
pd.list <- ColumnPlot(pettko.genes[, unique(gene.id)], return.plot.data = TRUE)
pb <- pd.list$pb

# merge family id
setkey(pd.list$plot.data, gene)
setkey(pettko.genes, gene.id)
plot.data <- pettko.genes[pd.list$plot.data]

# draw the plot
ggplot(plot.data, aes(y = symbol, x = log2FoldChange,
                      colour = log(mean.tpm + 0.5, 2))) +
  facet_grid(class ~ accession, scales = "free_y", space = "free_y") +
  ylab(NULL) +
  xlab(expression(L[2]*FC["PBM–SM"] %+-% "se ("*italic(n) == "3)")) +
  #    scale_colour_hue(c = 100, l = 50, h.start = 359) +
  scale_colour_gradientn(colours = heatscale,
                         limits = c(0, 10),
                         name = expression(Log[2]*TPM)) +
  geom_vline(xintercept = 0, size = 0.5, colour = "grey") +
  geom_vline(xintercept = c(-log(1.5, 2), log(1.5, 2)),
             size = 0.5, linetype = 2, colour = "grey") +
  geom_errorbarh(aes(xmax = log2FoldChange + lfcSE,
                     xmin = log2FoldChange - lfcSE),
                 height = 0.3, size = 0.3, colour = "black") +
  geom_point(size = 2)

################
# PLOT THE TPM #
################

# merge TPM values
pd2 <- tpm.long[pettko.genes]



ggplot(pd2[call == TRUE & class == "CYCD"],
       aes(x = accession, y = tpm, shape = accession, colour = stage,
                group = accession)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30)) +
  facet_wrap(~ gene, scales = "free") +
  geom_smooth(mapping = aes(group = stage), method = "loess", se = FALSE) +
  geom_point(position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1))





