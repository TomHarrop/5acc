#!/usr/bin/env Rscript

library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

# load lfc between stages per species
stage.results.file <- "output/deseq2/wald_stage_species/results_table.Rds"

if (!file.exists(stage.results.file)) {
  stop("Couldn't find results_table.Rds")
}
stage.results.table <- readRDS(stage.results.file)

# fix accession names
ord <- c(rufipogon = "O. rufipogon", indica = "O. sativa indica",
         japonica = "O. sativa japonica", barthii = "O. barthii",
         glaberrima = "O. glaberrima")
stage.results.table[ , accession := factor(plyr::revalue(
  accession,
  replace = ord),  levels = ord)]
setkey(stage.results.table, gene)

# select genes (testing)
# alog <- c('LOC_Os07g04670', 'LOC_Os02g07030', 'LOC_Os06g46030',
#           'LOC_Os02g41460', 'LOC_Os04g43580', 'LOC_Os10g33780',
#           'LOC_Os02g56610', 'LOC_Os01g61310', 'LOC_Os05g39500',
#           'LOC_Os05g28040')
# pd <- stage.results.table[alog]
tfdb.os <- readRDS("data/tfdb/tfdb_os.Rds")
goi <- tfdb.os[Family == "HB", unique(Protein.ID)]
tfdb.os[,unique(Family)]


goi <- unique(rownames(oryzr::SearchByGeneSymbol("RCN")))

pd <- stage.results.table[goi]

# add symbols
pd[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
pd[is.na(symbol), symbol := gene]

# cluster the y-axis
yclust.dt <- dcast(pd, symbol ~ accession, value.var = "log2FoldChange")
yclust <- as.matrix(data.frame(yclust.dt, row.names = "symbol"))
y.hc <- hclust(dist(yclust, method = "minkowski"), method = "ward.D2")
y.ord <- rev(rownames(yclust)[y.hc$order])
pd[, symbol := factor(symbol, levels = y.ord)]

ggplot(pd, aes(
  y = symbol, x = log2FoldChange, colour = symbol,
  xmax = log2FoldChange + lfcSE, xmin = log2FoldChange - lfcSE)) +
  theme_minimal() +
  facet_grid(~ accession) +
  scale_colour_hue(c = 100, l = 50, h.start = 359) +
  guides(colour = FALSE) +
  geom_vline(xintercept = 0, size = 0.5, colour = "grey") +
  geom_vline(xintercept = c(-log(1.5, 2), log(1.5, 2)),
             size = 0.5, linetype = 2, colour = "grey") +
  geom_errorbarh(height = 0.3, size = 0.3, colour = "black") +
  geom_point(size = 2) 

