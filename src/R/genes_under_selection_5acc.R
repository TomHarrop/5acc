#!/usr/bin/env Rscript

library(data.table)

# load data
gene.p.matrix <- readRDS("output/gene_p_matrix/gene_p_matrix.RDS")
annotated.window.genes <- readRDS(
  "output/genes_under_selection/annotated_window_genes.RDS")

p.value.columns <- names(gene.p.matrix)[
  !names(gene.p.matrix) %in% c("gene", "crowell_qtl", "lmi_qtl")]

# merge
gene.p.matrix[, symbols := NULL]
setkey(gene.p.matrix, gene)
setkey(annotated.window.genes, gene)
selection.with.p <- gene.p.matrix[annotated.window.genes]

setcolorder(selection.with.p, c(
  "Chr", "window.start", "window.end", "test", "gene", "symbol",
  p.value.columns, "crowell_qtl", "lmi_qtl",
  "annotation", "OgroObjective", "OgroRef"))

# extract significant genes
any.sig <- selection.with.p[!is.na(domestication),
                            .(any.sig = any(.SD < 0.1)),
                            by = gene,
                            .SDcols = p.value.columns][, unique(gene)]
sig.in.selection <- selection.with.p[gene %in% any.sig |
                                       crowell_qtl == TRUE |
                                       lmi_qtl == TRUE]

# write output
outdir = "output/genes_under_selection"
if(!dir.exists(outdir)){
  dir.create(outdir)
}
saveRDS(selection.with.p, paste0(outdir, "/window_genes_x_5acc.RDS"))
write.table(selection.with.p,
            paste0(outdir, "/window_genes_x_5acc_all.tab"),
            sep = "\t",
            quote = FALSE,
            na = "",
            row.names = FALSE)
write.table(sig.in.selection,
            paste0(outdir, "/window_genes_x_5acc_sig.tab"),
            sep = "\t",
            quote = FALSE,
            na = "",
            row.names = FALSE)

sinfo <- rutils::GitSessionInfo()
writeLines(sinfo, paste0(outdir, "/SessionInfo.genes_under_selection_5acc.txt"))





