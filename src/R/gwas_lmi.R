#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(gtable)
library(gridExtra)

rutils::GenerateMessage("Annotate LMI GWAS regions")

# load l2fc data
rutils::GenerateMessage("Loading DE results")
stage.l2fc.dom.padj <- readRDS(
  "output/deseq2/wald_tests/stage_l2fc_dom_padj.Rds")

# load tpm data
rutils::GenerateMessage("Load TPM gene expression data")
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

# load gwas qtl regions
rutils::GenerateMessage("Load LMI GWAS results")
lmi.gwas <- data.table(read.csv("data/gwas/lmi/QTL_VN_coordinates.csv"))

# setup bin.name
rutils::GenerateMessage("Munge regions")
setkey(lmi.gwas, Chr, Start, End)
lmi.gwas[, bin.name := paste0(
  paste(unique(Trait), collapse = "+"), "|",
  paste(unique(Place), collapse = "+"), "|",
  paste(unique(Group), collapse = "+"), "|",
  "Chr", Chr, ":", Start, "-", End), by = .(Chr, Start, End)]
lmi.gwas[, bin.name := factor(bin.name, levels = unique(bin.name))]
setkey(lmi.gwas, bin.name)
lmi.gwas <- unique(lmi.gwas)

# get genes from phytozome:
# IRD network blocks GET after phytozome returns 302 from initial POST. Tell
# RCurl to use POST instead.
rutils::GenerateMessage("Setting up biomaRt")
options(RCurlOptions = list(followlocation = TRUE, postredir = 2L))

# set up biomaRt
phytozome <- biomaRt::useMart(biomart = "phytozome_mart",
                              dataset = "phytozome",
                              host = "phytozome.jgi.doe.gov",
                              path = "/biomart/martservice/")

# setup query
bm.filters <- c("organism_id", "region_name", "region_start", "region_end")
bm.attributes <- c("gene_name1")

# setup query function
curl.handle <- RCurl::getCurlHandle()
curl.handle <- RCurl::curlSetOpt(
  .opts = list(followlocation = TRUE, postredir = 2L), curl = curl.handle)
RetrieveRegionGenes <- function(region_name, region_start, region_end) {
  rutils::PrintF("biomaRt query: %s:%s-%s\n",
                 region_name, region_start, region_end)
  bm.values <- list(
    organism_id = "323", # O. sativa
    region_name = region_name,
    region_start = region_start,
    region_end = region_end)
  unique(biomaRt::getBM(bm.attributes, bm.filters, bm.values,
                        phytozome, uniqueRows = TRUE,
                        curl = curl.handle)$gene_name1)
}

# run query
rutils::GenerateMessage("Running biomaRt queries per QTL")
lmi.gwas[, region.name := paste0("Chr", Chr)]
qtl.to.gene <- lmi.gwas[, .(
  gene = RetrieveRegionGenes(region.name, Start, End)),
  by = bin.name]

# weed out strange answers
rutils::GenerateMessage("Merging results")
qtl.to.gene <- qtl.to.gene[!is.na(gene)][
  grep("LOC_Os[[:digit:]]+g[[:digit:]]+", gene)]

# merge results
setkey(lmi.gwas, bin.name)
setkey(qtl.to.gene, bin.name)
gwas.lmi.genes <- qtl.to.gene[lmi.gwas, .(
  bin.name, gene, Chr, region.name, Start, End, Trait, Place, Group
)]

# collapse gwas.lmi.genes results by gene
setkey(gwas.lmi.genes, Chr, Start, End, Place)
gwas.lmi.genes[, gene.bin := paste(bin.name, collapse = "\n"), by = gene]
gwas.lmi.genes[, gene.bin := factor(gene.bin, levels = unique(gene.bin))]

# join domestication results to gwas
rutils::GenerateMessage("Joining domestication results")
setkey(stage.l2fc.dom.padj, gene, accession)
setkey(gwas.lmi.genes, gene)
gwas.with.p <- stage.l2fc.dom.padj[gwas.lmi.genes]

# tidy (select rows)
setkey(gwas.with.p, gene, accession)
gwas.with.p <- unique(gwas.with.p)[, .(
  gene, dom_all, dom_africa, dom_asia, dom_japonica, dom_indica, accession,
  baseMean, log2FoldChange, lfcSE, padj, gene.bin, region.name, Start, End)]

# annotate and output
rutils::GenerateMessage("Adding annotations")
gwas.annot <- gwas.with.p[!is.na(gene), data.table(oryzr::LocToGeneName(gene),
                                                   keep.rownames = TRUE)]
setnames(gwas.annot, "rn", "gene")
setkey(gwas.annot, gene)
gwas.annot <- unique(gwas.annot)
gwas.annot.full <- gwas.annot[gwas.with.p]
setkey(gwas.annot.full, gene.bin, gene, accession)
my.names <- c("region.name", "Start", "End", "gene.bin", "gene", "symbols",
              "names", "MsuAnnotation", "dom_all", "dom_africa", "dom_asia",
              "dom_japonica", "dom_indica", "accession", "baseMean",
              "log2FoldChange", "lfcSE", "padj", "OgroObjective", "OgroRef",
              "RapID")
setcolorder(gwas.annot.full, my.names)
gwas.dom.sig <- gwas.annot.full[
  dom_all < 0.05 | dom_africa < 0.05 | dom_asia < 0.05 | dom_japonica < 0.05 |
    dom_indica < 0.05
  ]

# plot the regions
PlotGwasBin <- function(bin) {
  
  rutils::PrintF("gene.bin: %s\n", bin)
  
  # subset plot data
  plot.data <- gwas.annot.full[gene.bin == bin]
  plot.title <- plot.data[, unique(gene.bin)]
  
  # if no genes, return empty plot
  if (dim(plot.data[!is.na(accession)])[1] == 0) {
    return(ggplot() + ggtitle(plot.title))
  }
  
  # clean up NAs
  plot.data <- plot.data[!is.na(accession)]
  
  # add tpm data
  setkey(tpm.basemean, gene, accession)
  setkey(plot.data, gene, accession)
  plot.data <- tpm.basemean[plot.data]
  
  # fix gene names
  plot.data[is.na(symbols), symbols := gene]
  
  # sort y-axis by lfc in rufipogon
  setkey(plot.data, accession, log2FoldChange)
  y.ord <- plot.data[accession == levels(accession)[1],
                     as.character(unique(symbols))]
  id.ord <- plot.data[accession == levels(accession)[1],
                      as.character(unique(gene))]
  plot.data[, symbols := factor(symbols, levels = rev(y.ord))]
  
  # colour the points
  heatscale <- RColorBrewer::brewer.pal(6, "YlOrRd")
  
  # main plot
  main.plot <- ggplot(plot.data, aes(
    y = symbols, x = log2FoldChange, colour = mean.tpm + 1,
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
  
  # convert p < 0.05 domestication results to numerical for symbols
  tt.gene.dom <- plot.data[,lapply(.SD, function(x) as.character(x <  0.05)),
                           .SDcols = c(
                             "dom_all", "dom_africa", "dom_asia", "dom_indica",
                             "dom_japonica"),
                           by = gene]
  tt.gene.dom[dom_all == "TRUE", dom_all := "All"]
  tt.gene.dom[dom_asia == "TRUE", dom_asia := "Asia"]
  tt.gene.dom[dom_africa == "TRUE", dom_africa := "Africa"]
  tt.gene.dom[dom_indica == "TRUE", dom_indica := "Indica"]
  tt.gene.dom[dom_japonica == "TRUE", dom_japonica := "Japonica"]
  tt.long <- melt(tt.gene.dom, id.vars = "gene", variable.name = "domestication",
                  value.name = "shape")
  tt.long[shape == "FALSE", shape := NA]
  shape.levels <- c("All", "Africa", "Asia", "Indica",
    "Japonica")
  tt.long[, shape := factor(shape, levels = shape.levels)]
  
  # if no significant genes, just return main plot
  if (dim(tt.long[!is.na(shape)])[1] == 0) {
    return(ggplotGrob(main.plot))
  }
  
  # subset by plotted genes
  setkey(plot.data, gene)
  setkey(tt.long, gene)
  symbol.plot.data <- tt.long[plot.data, .(
    gene, domestication, shape, symbols
  ), by = .EACHI]
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
      breaks = shape.levels,
      labels = shape.levels,
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

  return(plot.complete)

}

rutils::GenerateMessage("Generating plots")
# one plot per bin in a list
all.bins <- gwas.annot.full[!is.na(gene), unique(as.character(gene.bin))]
list.of.plots <- lapply(all.bins, PlotGwasBin)

# print to pdf
plots.arranged <- do.call(gridExtra::marrangeGrob, args = list(
  grobs = list.of.plots,
  nrow = 1, ncol = 1))

# save output
out.dir <- "output/gwas/lmi"
rutils::GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

write.table(gwas.annot.full[, setdiff(my.names, "gene.bin"), with = FALSE],
            paste0(out.dir, "/gwas_lmi_all.tab"),
            quote = FALSE, sep = "\t",
            na = "", row.names = FALSE)
write.table(gwas.dom.sig[, setdiff(my.names, "gene.bin"), with = FALSE],
            paste0(out.dir, "/gwas_lmi_sig.tab"),
            quote = FALSE, sep = "\t",
            na = "", row.names = FALSE)

saveRDS(gwas.annot.full, paste0(out.dir, "/gwas_annot_full.Rds"))
saveRDS(gwas.lmi.genes, paste0(out.dir, "/gwas_lmi_genes.Rds"))

rutils::GenerateMessage("Writing plots")
ggsave(paste0(out.dir, "/gwas_lmi_rnaseq.pdf"), plots.arranged,
       width = 10.69, height = 7.27)

# save logs
sInf <- rutils::GitSessionInfo()
logLocation <- paste0(out.dir, "/SessionInfo.gwas_lmi.txt")
writeLines(sInf, logLocation)

rutils::GenerateMessage("Done")

quit(save = "no", status = 0)
