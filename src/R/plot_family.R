#!/usr/bin/env Rscript

library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

# load lfc between stages per species
stage.results.file <-
  "output/deseq2/wald_tests/stage_species_results_table.Rds"
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
alog <- c('LOC_Os07g04670', 'LOC_Os02g07030', 'LOC_Os06g46030',
          'LOC_Os02g41460', 'LOC_Os04g43580', 'LOC_Os10g33780',
          'LOC_Os02g56610', 'LOC_Os01g61310', 'LOC_Os05g39500',
          'LOC_Os05g28040')
pd <- stage.results.table[alog]

goi <- c('LOC_Os01g44170', 'LOC_Os01g49310', 'LOC_Os01g65370', 'LOC_Os01g67310', 'LOC_Os02g13800', 'LOC_Os02g16690', 'LOC_Os03g42464', 'LOC_Os04g30610', 'LOC_Os04g40630', 'LOC_Os04g55159', 'LOC_Os04g57550', 'LOC_Os04g59260', 'LOC_Os05g08420', 'LOC_Os05g15040', 'LOC_Os05g29810', 'LOC_Os05g37520', 'LOC_Os05g41760', 'LOC_Os05g41780', 'LOC_Os06g24910', 'LOC_Os06g42560', 'LOC_Os08g19670', 'LOC_Os08g31580', 'LOC_Os08g34640', 'LOC_Os08g39450', 'LOC_Os09g36420', 'LOC_Os10g38489', 'LOC_Os11g08440', 'LOC_Os12g31748', 'LOC_Os12g32660')




tfdb.os <- readRDS("data/tfdb/tfdb_os.Rds")
tfdb.os[,unique(Family)]

# goi <- unique(rownames(oryzr::SearchByGeneSymbol("RCN")))
LOB
SBP
GRF
GRAS
ARF
AUX/IAA

goi <- tfdb.os[Family == "SBP", unique(Protein.ID)]

# select data for plot
pd <- stage.results.table[goi]
pd <- pd[!is.na(accession)]

# add symbols
pd[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
pd[is.na(symbol), symbol := gene]

# cluster the y-axis, not really nice
# yclust.dt <- dcast(pd, symbol ~ accession, value.var = "log2FoldChange")
# yclust <- as.matrix(data.frame(yclust.dt, row.names = "symbol"))
# y.hc <- hclust(dist(yclust, method = "minkowski"), method = "ward.D2")
# y.ord <- rev(rownames(yclust)[y.hc$order])
# pd[, symbol := factor(symbol, levels = y.ord)]

# sort y-axis by lfc in rufipogon
setkey(pd, "accession", "log2FoldChange")
y.ord <- pd[accession == levels(accession)[1], as.character(unique(symbol))]
pd[, symbol := factor(symbol, levels = rev(y.ord))]

family_plot <- ggplot(pd, aes(
  y = symbol, x = log2FoldChange, colour = symbol,
  xmax = log2FoldChange + lfcSE, xmin = log2FoldChange - lfcSE)) +
  theme_grey(base_size = 16) +
  theme(axis.text.y	= element_text(face = "italic"),
        strip.text = element_text(face = "italic")) +
  xlab(expression(L[2]*FC["(PBM–SM)"] %+-% se)) + ylab(NULL) +
  facet_grid(~ accession) +
  scale_colour_hue(c = 100, l = 50, h.start = 359) +
  guides(colour = FALSE) +
  geom_vline(xintercept = 0, size = 0.5, colour = "grey") +
  geom_vline(xintercept = c(-log(1.5, 2), log(1.5, 2)),
             size = 0.5, linetype = 2, colour = "grey") +
  geom_errorbarh(height = 0.1, size = 0.5, colour = "black") +
  geom_point(size = 2) 

cairo_pdf("~/Desktop/alog.pdf", width = 10, height = 7.5, pointsize = 16,
          bg = "transparent")
print(family_plot)
dev.off()