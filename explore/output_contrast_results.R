#!/usr/bin/env Rscript

library(data.table)

# function to annotate the tables
AnnotateResults <- function(x) {
  my.dt <- copy(x)
  setkey(my.dt, "gene")
  my.dt[, c(
    'RapID', 'symbols', 'names', 'MsuAnnotation', 'OgroObjective',
    'OgroRef') :=
      oryzr::LocToGeneName(gene),
    by = gene]
  return(my.dt)
}

# 1. species differential expression
species.results <- readRDS("output/deseq2/wald_species/contrast_results.Rds")
species.results[, unique(contrast)]
sig.species.results <- species.results[padj < 0.05 & abs(log2FoldChange) > 1]
sig.species.results <- AnnotateResults(sig.species.results)
setorder(sig.species.results, contrast, -log2FoldChange)

# 2. stage differential expression
# 2a. all species
stage.results <- readRDS("output/deseq2/wald_species/stage_results_table.Rds")
sig.stage.results <- stage.results[padj < 0.05 & abs(log2FoldChange) > log(1.5, 2)]
sig.stage.results <- AnnotateResults(sig.stage.results)
sig.stage.results[, abs.l2fc := abs(log2FoldChange)]
setorder(sig.stage.results, -abs.l2fc)
sig.stage.results[, abs.l2fc := NULL]
# 2b. within continents
stage.continent.results <- readRDS(
  "output/deseq2/wald_stage_continent/results_table.Rds")
sig.sc.results <- stage.continent.results[
  padj < 0.05 & abs(log2FoldChange) > log(1.5, 2)]
sig.sc.results <- AnnotateResults(sig.sc.results)
sig.sc.results[, abs.l2fc := abs(log2FoldChange)]
setorder(sig.sc.results, continent, -abs.l2fc)
sig.sc.results[, abs.l2fc := NULL]

# 3. stage:domestication interaction
sd.results <- readRDS("output/deseq2/wald_domestication/results_table.Rds")
sig.sd.results <- sd.results[padj < 0.05]
sig.sd.results <- AnnotateResults(sig.sd.results)
sig.sd.results[, abs.l2fc := abs(log2FoldChange)]
setorder(sig.sd.results, -abs.l2fc)
sig.sd.results[, abs.l2fc := NULL]

# write output
WriteOutput <- function(x, file.name, out.dir = "explore/xlsx") {
  write.table(x, paste0(out.dir, "/", file.name), quote = FALSE, sep = "\t",
              na = "", row.names = FALSE)
}

WriteOutput(sig.species.results, "species_diff_exp.tsv")
WriteOutput(sig.stage.results, "stage_diff_exp.tsv")
WriteOutput(sig.sd.results, "stage_domestication_interaction.tsv")
WriteOutput(sig.sc.results, "stage_by_continent.tsv")

# tpm <- readRDS("output/tpm/tpm_with_calls.Rds")
# 
# 
# ggplot(tpm["LOC_Os03g59550"],
#        aes(x = stage, y = tpm, colour = species, group = species)) +
#   scale_colour_brewer(palette = "Set1") +
#   geom_smooth(method = "lm", se = FALSE) +
#   geom_point(alpha = 0.75, position = position_jitter(width = 0.1))
