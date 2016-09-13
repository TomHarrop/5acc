library(data.table)
library(ggplot2)
library(gtable)
library(gridExtra)

# load gwas / domestication table
gwas.dom <- readRDS("output/gwas/gwas_domestication.Rds")

# load tpm data
tpm.long <- readRDS("output/tpm/tpm.Rds")
ord <- c("O. rufipogon", "O. sativa indica", "O. sativa japonica",
         "O. barthii", "O. glaberrima")
tpm.long[, accession := substr(sample, 1, 1)]
tpm.long[ , accession := factor(plyr::mapvalues(
  accession,
  from = c("R", "I", "J", "B", "G"),
  to = ord),  levels = ord)
  ]
tpm.basemean <- tpm.long[, .(mean.tpm = mean(tpm)), by = .(gene, accession)]

# make numeric bin id
gwas.dom[, chromosome.number := as.numeric(
  gsub("Chr", "", region_name, fixed = TRUE))]
setkey(gwas.dom, chromosome.number, region_start, region_end)
gwas.dom[, bid.num := as.numeric(
  factor(Bin_id, levels = as.character(unique(Bin_id))))]

# test plot with some data subset
gwas.dom[dom_all < 0.05 | dom_africa < 0.05 | dom_asia < 0.05 |
           dom_indica < 0.05 | dom_japonica < 0.05, length(unique(gene)),
         by = bid.num]
gwas.dom[dom_asia < 0.05, length(unique(gene)), by = bid.num]
gwas.dom[Bin_id == "chr4: 30880000-31310000"] # chr 4 superlocus fig 7
gwas.dom[, unique(Bin_id)]

# check:LOC_Os02g47660, bid 41

PlotGwasBin <- function(bid.number) {
  
  cat("bid.number:", bid.number, "\n")
  
  plot.data <- gwas.dom[bid.num == bid.number]

  # plot.title
  plot.title <- plot.data[, paste(
    unique(Bin_id),
    unique(trait.collapsed),
    unique(subpop.collapsed),
    sep = " | "
  )]
  
  if (dim(plot.data[!is.na(accession)])[1] == 0) {
    return(ggplot() + ggtitle(plot.title))
  }
  
  # add tpm data
  setkey(tpm.basemean, gene, accession)
  setkey(plot.data, gene, accession)
  plot.data <- tpm.basemean[plot.data]
  
  # add symbols
  plot.data <- plot.data[!is.na(accession)]
  plot.data[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
  plot.data[is.na(symbol), symbol := gene]
  
  # sort y-axis by lfc in rufipogon
  setkey(plot.data, accession, log2FoldChange)
  y.ord <- plot.data[accession == levels(accession)[1],
                     as.character(unique(symbol))]
  id.ord <- plot.data[accession == levels(accession)[1],
                      as.character(unique(gene))]
  plot.data[, symbol := factor(symbol, levels = rev(y.ord))]
  
  # colour the points
  heatscale <- RColorBrewer::brewer.pal(6, "YlOrRd")
  
  # main plot
  main.plot <- ggplot(plot.data, aes(
    y = symbol, x = log2FoldChange, colour = mean.tpm + 1,
    xmax = log2FoldChange + lfcSE, xmin = log2FoldChange - lfcSE)) +
    theme_grey(base_size = 12) +
    theme(axis.text.y	= element_text(face = "italic"),
          strip.text.x = element_text(face = "italic"),
          legend.title = element_text(size = rel(0.5)),
          legend.text = element_text(size = rel(0.5)),
          legend.key.size = unit(0.8, "lines"),
          plot.margin = unit(c(1, 1, 1, 1), "lines")) +
    xlab(expression(L[2]*FC[(PBM -  SM)] %+-% se)) + ylab(NULL) +
    facet_grid(~ accession, scales = "free_y", space = "free_y") +
    scale_colour_gradientn(
      limits = c(1, 2^10),
      colours = heatscale,
      trans = scales::log2_trans(),
      breaks = scales::trans_breaks("log2", function(x) 2^x),
      name = expression(bar(TPM) + 1),
      na.value = "#800026") +
    geom_vline(xintercept = 0, size = 0.5, colour = "grey") +
    geom_vline(xintercept = c(-log(1.5, 2), log(1.5, 2)),
               size = 0.5, linetype = 2, colour = "grey") +
    geom_errorbarh(height = 0.1, size = 0.5, colour = "black") +
    geom_point(size = 2) +
    ggtitle(plot.title)
  
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
    theme_void(base_size = 12) +
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
  
  # add the panel to main.figure
  main.plot.grob <- ggplotGrob(main.plot)
  axis.l <- subset(main.plot.grob$layout, name == "axis-l")
  plot.with.panel <- gtable_add_cols(
    main.plot.grob,
    unit(2/3, "null"),
    axis.l$l - 1)
  plot.with.panel <- gtable_add_grob(
    plot.with.panel, symbol.panel,
    t = axis.l$t, b = axis.l$b,
    l = axis.l$l, r = axis.l$r)
  
  # add the legend outside the plot
  plot.complete <- gtable_add_cols(
    plot.with.panel,
    unit(2/3, "null"),
    axis.l$r - 1
  )
  plot.complete <- gtable_add_grob(
    plot.complete, symbol.legend,
    t = axis.l$t, b = axis.l$b,
    l = axis.l$r, r = axis.l$r,
    clip = "on"
  )
  
  # gtable_show_layout(plot.complete)
  # gtable_show_layout(plot.with.panel)
  
  return(plot.complete)
}

# get plots, one per bin
all.qtls <- gwas.dom[, unique(bid.num)]
all.qtls <- all.qtls[order(all.qtls)]
list.of.plots <- lapply(all.qtls, PlotGwasBin)

# print to pdf
plots.arranged <- do.call(gridExtra::marrangeGrob, args = list(
  grobs = list.of.plots,
  nrow = 1, ncol = 1))

ggsave("explore/gwas_dom.pdf", plots.arranged, width = 10.69, height = 7.27)
