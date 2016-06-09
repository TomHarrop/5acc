#!/usr/bin/env Rscript

library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}
GenerateMessage("Look for domestication-associated genes in GWAS regions")

# load gwas data
GenerateMessage("Loading genes from GWAS regions")
gwas.crowell.genes.file <- "output/gwas/crowell/gwas_crowell_genes.Rds"
if (!file.exists(gwas.crowell.genes.file)) {
  stop("Couldn't find gwas_crowell_genes.Rds")
}
gwas.crowell.genes <- readRDS(gwas.crowell.genes.file)

# load domestication results
# TODO: TIDY UP LOADING, DO MUNGING IN extract_gene_lists.R
GenerateMessage("Loading domestication Wald test results")
dom.asia <- readRDS(
  "output/deseq2/wald_tests/domestication_asia_results_table.Rds")
dom.continent <- readRDS(
  "output/deseq2/wald_tests/domestication_by_continent_results_table.Rds")
dom.all <- rbind(dom.asia, dom.continent)
p.by.dom <- dcast(dom.all, gene ~ domestication, value.var = "padj")

# join domestication to gwas
setkey(p.by.dom, gene)
setkey(gwas.crowell.genes, gene)
gwas.with.p <- p.by.dom[gwas.crowell.genes]

# set up a plot (do with plot_family.R)
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

# select data for plot
# goi <- gwas.with.p[asia < 0.05 | africa < 0.05 | indica < 0.05 |
#                      japonica < 0.05,.(
#                        gene, africa, asia, indica, japonica, Bin_id, TRAIT,
#                        SUB_POP
#                      )]
goi <- gwas.with.p[asia < 0.05 | indica < 0.05 |
                     japonica < 0.05,.(
                       gene, africa, asia, indica, japonica, Bin_id, TRAIT,
                       SUB_POP
                     )]



# collapse genes by TRAIT and SUB_POP
goi[, trait.collapsed := paste(unique(TRAIT), collapse = ", "), by = gene]
goi[, subpop.collapsed := paste(unique(SUB_POP), collapse = ", "), by = gene]
setkey(goi, gene)
goi <- unique(goi[, .(gene, africa, asia, indica, japonica, Bin_id,
                      trait.collapsed, subpop.collapsed)])
goi[, facet.y := paste(Bin_id, trait.collapsed, subpop.collapsed)]

# make plot
pd <- stage.results.table[goi]
pd <- pd[!is.na(accession)]

# add symbols
pd[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
pd[is.na(symbol), symbol := gene]

# sort y-axis by lfc in rufipogon
setkey(pd, "accession", "log2FoldChange")
y.ord <- pd[accession == levels(accession)[1], as.character(unique(symbol))]
pd[, symbol := factor(symbol, levels = rev(y.ord))]

#pd <- pd[symbol %in% c("ERF77", "SPL10")]
# diff bin ID?
pd[, bid.num := as.numeric(factor(Bin_id))]
pd2 <- pd[bid.num < 31 & bid.num > 20]

ggplot(pd2, aes(x = log2FoldChange, y = symbol)) +
  facet_grid(bid.num ~ accession, scales = "free_y", space = "free_y") +
  geom_point()


ggplot(pd2, aes(
  y = symbol, x = log2FoldChange, colour = symbol,
  xmax = log2FoldChange + lfcSE, xmin = log2FoldChange - lfcSE)) +
  theme_grey(base_size = 16) +
  theme(axis.text.y	= element_text(face = "italic"),
        strip.text = element_text(face = "italic")) +
  xlab(expression(L[2]*FC["(PBM–SM)"] %+-% se)) + ylab(NULL) +
  facet_grid(bid.num ~ accession, scales = "free_y", space = "free_y") +
  scale_colour_hue(c = 100, l = 50, h.start = 359) +
  guides(colour = FALSE) +
  geom_vline(xintercept = 0, size = 0.5, colour = "grey") +
  geom_vline(xintercept = c(-log(1.5, 2), log(1.5, 2)),
             size = 0.5, linetype = 2, colour = "grey") +
  geom_errorbarh(height = 0.1, size = 0.5, colour = "black") +
  geom_point(size = 2) 
