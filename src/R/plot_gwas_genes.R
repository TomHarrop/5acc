library(data.table)
library(ggplot2)
library(gtable)

gwas.dom <- readRDS("output/gwas/gwas_domestication.Rds")

# make numeric bin id
gwas.dom[, bid.num := as.numeric(factor(Bin_id))]

# test plot with some data subset
gwas.dom[dom_all < 0.05 | dom_africa < 0.05 | dom_asia < 0.05 |
            dom_indica < 0.05 | dom_japonica < 0.05, length(unique(gene)),
         by = bid.num]
gwas.dom[dom_asia < 0.05, length(unique(gene)), by = bid.num]
gwas.dom[Bin_id == "chr4: 30880000-31310000"] # chr 4 superlocus fig 7
gwas.dom[, unique(Bin_id)]
plot.data <- gwas.dom[bid.num == 121]
plot.data <- plot.data[dom_all < 0.05 | dom_africa < 0.05 | dom_asia < 0.05 |
                         dom_indica < 0.05 | dom_japonica < 0.05]

# add symbols
plot.data[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
plot.data[is.na(symbol), symbol := gene]

# sort y-axis by lfc in rufipogon
setkey(plot.data, "accession", "log2FoldChange")
y.ord <- plot.data[accession == levels(accession)[1],
                   as.character(unique(symbol))]
id.ord <- plot.data[accession == levels(accession)[1],
                as.character(unique(gene))]
plot.data[, symbol := factor(symbol, levels = rev(y.ord))]

# colour the points
plot.cols <- plot.data[, wesanderson::wes_palette(
  "Zissou", length(unique(symbol)), "continuous")]

# main plot
main.plot <- ggplot(plot.data, aes(
  y = symbol, x = log2FoldChange, colour = symbol,
  xmax = log2FoldChange + lfcSE, xmin = log2FoldChange - lfcSE)) +
  theme_grey(base_size = 16) +
  theme(axis.text.y	= element_text(face = "italic"),
        strip.text = element_text(face = "italic")) +
  xlab(expression(L[2]*FC["(PBM–SM)"] %+-% se)) + ylab(NULL) +
  facet_grid(bid.num ~ accession, scales = "free_y", space = "free_y") +
  #scale_colour_hue(c = 100, l = 50, h.start = 359) +
  scale_colour_manual(values = plot.cols) +
  guides(colour = FALSE) +
  geom_vline(xintercept = 0, size = 0.5, colour = "grey") +
  geom_vline(xintercept = c(-log(1.5, 2), log(1.5, 2)),
             size = 0.5, linetype = 2, colour = "grey") +
  geom_errorbarh(height = 0.1, size = 0.5, colour = "black") +
  geom_point(size = 2) 

# make domestication key

# convert p < 0.05 domestication results to numerical for symbols
tt.gene.dom <- gwas.dom[,lapply(.SD, function(x) as.numeric(x <  0.05)),
                        .SDcols = c(
                          "dom_all", "dom_africa", "dom_asia", "dom_indica",
                          "dom_japonica"),
                        by = gene]
tt.gene.dom[dom_africa == 1, dom_africa := 2]
tt.gene.dom[dom_asia == 1, dom_asia := 3]
tt.gene.dom[dom_indica == 1, dom_indica := 4]
tt.gene.dom[dom_japonica == 1, dom_japonica := 5]
tt.long <- melt(tt.gene.dom, id.vars = "gene", variable.name = "domestication",
                value.name = "shape")
tt.long[shape == 0, shape := NA]
tt.long[, shape := factor(shape)]

# subset by plotted genes
setkey(plot.data, gene, accession)
setkey(tt.long, gene)
symbol.plot.data <- tt.long[plot.data, .(
  gene, domestication, shape, symbol
)]
setkey(symbol.plot.data, gene, domestication)
symbol.plot.data <- unique(symbol.plot.data)

# order tt.gene by plot.data symbols
symbol.plot.data[, gene := factor(gene, levels = rev(id.ord))]

# draw the symbol panel
symbol.plot <- ggplot(symbol.plot.data,
                       aes(y = gene, x = domestication, shape = shape)) +
  theme_void() +
#  theme(legend.title = element_text()) +
  geom_point(size = 1, colour = "red") +
  scale_shape_manual(
    name = "Domestication",
    breaks = symbol.plot.data[,levels(shape)],
    labels = c("All", "Africa", "Asia", "Indica",
               "Japonica"),
    values = c(0, 2, 3, 4, 5),
    drop = FALSE, na.value = NA)

# extract panel from symbol.panel
symbol.grob <- ggplotGrob(symbol.plot)
symbol.panel <- gtable_filter(symbol.grob, "panel")
symbol.legend <- gtable_filter(symbol.grob, "guide-box")

# add the panel main.figure
main.plot.grob <- ggplotGrob(main.plot)
strip.right <- subset(main.plot.grob$layout, name == "strip-right")
strip.right <- subset(main.plot.grob$layout, name == "axis-l")
plot.with.panel <- gtable_add_cols(
  main.plot.grob,
  unit(0.5, "null"),
  strip.right$l - 1)
plot.with.panel <- gtable_add_grob(
  plot.with.panel, symbol.panel,
  t = strip.right$t, b = strip.right$b,
  l = strip.right$l, r = strip.right$r)

# add the legend outside the plot
gtable_show_layout(plot.complete)
gtable_show_layout(plot.with.panel)

plot.complete <- gtable_add_cols(
  plot.with.panel,
  unit(2/3, "null"),
  strip.right$r - 1
)
plot.complete <- gtable_add_grob(
  plot.complete, symbol.legend,
  t = strip.right$t, b = strip.right$b,
  l = strip.right$r, r = strip.right$r,
  clip = "on"
)

grid::grid.newpage()
grid::grid.draw(plot.complete)
grid::grid.draw(plot.with.panel)




