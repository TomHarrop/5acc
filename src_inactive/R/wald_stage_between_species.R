#!/usr/bin/env Rscript

library(DESeq2)
library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

GenerateMessage("Differential expression within stages between species")

# set up parallel processing
cpus <- as.numeric(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
if (is.na(cpus)) {
  cpus <- 1
}
GenerateMessage(paste("Allocating", cpus, "cpu(s)"))
BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# load data
GenerateMessage("Loading base DESeq2 object")
dds.file <- "output/deseq2/dds.Rds"
if (!file.exists(dds.file)) {
  stop("Couldn't find dds.Rds")
}
dds <- readRDS(dds.file)

# pre-filter based on tpm cutoffs
GenerateMessage("Removing undetected genes")
detected.genes.file <- "output/tpm/detected_genes.Rds"
if (!file.exists(detected.genes.file)) {
  stop("Couldn't fine detected_genes.Rds")
}
detected.genes <- readRDS(detected.genes.file)

# re-run DESeq2 per group
col.data <- SummarizedExperiment::colData(dds)
col.data$group <- paste(col.data$accession, col.data$stage, sep = "_")
dds_group <- DESeqDataSetFromMatrix(counts(dds), col.data, ~group)
dds_group <- DESeq(dds_group, parallel = TRUE)

# set up contrasts
all_species_contrasts <- data.table(
  expand.grid(s1 = levels(col.data$accession),
              s2 = levels(col.data$accession))
)
all_species_contrasts[, s1 := as.character(s1)]
all_species_contrasts[, s2 := as.character(s2)]
spec_contrasts <- unique(all_species_contrasts[
  , .(comb = paste(sort(c(s1, s2)), collapse = "_")),
  by = .(s1, s2)],
  by = "comb")[, .(s1, s2)][s1 != s2]

# function for getting contrast results
ExtractContrastResults <- function(denom,
                                   num,
                                   deseq_object,
                                   stage,
                                   alpha,
                                   lfcThreshold) {
  my_contrast <- c(
    "group",
    paste(num, stage, sep = "_"),
    paste(denom, stage, sep = "_")
  )
  GenerateMessage("Extracting results for contrast")
  message(paste0(my_contrast[[1]], ": ",
                 my_contrast[[2]], " vs. ",
                 my_contrast[[3]]))
  my_res <- results(object = deseq_object,
                    contrast = my_contrast,
                    lfcThreshold = lfcThreshold,
                    alpha = alpha,
                    parallel = TRUE)
  my_dt <- data.table(data.frame(my_res), keep.rownames = TRUE)
  setnames(my_dt, "rn", "gene")
  return(my_dt)
}


# run contrast results
contrast_alpha <- 0.05
contrast_lfc <- 1
contrast_results_list <- list(
  PBM = spec_contrasts[, ExtractContrastResults(num = s1,
                                                denom = s2,
                                                deseq_object = dds_group,
                                                stage = "PBM",
                                                alpha = contrast_alpha,
                                                lfcThreshold = contrast_lfc),
                       by = .(s1, s2)],
  SM = spec_contrasts[, ExtractContrastResults(num = s1,
                                               denom = s2,
                                               deseq_object = dds_group,
                                               stage = "SM",
                                               alpha = contrast_alpha,
                                               lfcThreshold = contrast_lfc),
                      by = .(s1, s2)])
contrast_results <- rbindlist(contrast_results_list, idcol = "stage")
setnames(contrast_results,
         c("s1", "s2"),
         c("wald_numerator", "wald_denominator"))
summary_table <- contrast_results[
  padj < contrast_alpha,
  .(alpha = contrast_alpha,
    lfcThreshold = contrast_lfc,
    number_of_genes = length(unique(gene))),
  by = .(stage, wald_numerator, wald_denominator)]

# save output
GenerateMessage("Saving output")
out.dir <- "output/deseq2/wald_tests"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(contrast_results, paste0(out.dir, "/stage_between_species_results_table.Rds"))
saveRDS(summary_table, paste0(out.dir, "/stage_between_species_summary_table.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.wald_stage_between_species.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)

# dirty venn diagram
ExtractContrastGenes <- function(s1, s2, stage) {
  contrast_results[padj < contrast_alpha &
                     (wald_numerator %in% c(s1, s2)) &
                     (wald_denominator %in% c(s1, s2)) &
                     stage == stage, unique(gene)]
}
pbm_bj <- ExtractContrastGenes("barthii", "japonica", "PBM")

contrast_results[padj < contrast_alpha &
                   (wald_numerator %in% c("barthii", "japonica")) &
                   (wald_denominator %in% c("barthii", "japonica")) &
                   stage == "PBM", unique(gene)]

Set1 <- RColorBrewer::brewer.pal(9, "Set1")
vd <- VennDiagram::venn.diagram(
  x = list(
    barthii = ExtractContrastGenes("barthii", "japonica", "PBM"),
    rufipogon = ExtractContrastGenes("rufipogon", "japonica", "PBM"),
    indica = ExtractContrastGenes("indica", "japonica", "PBM"),
    glaberrima = ExtractContrastGenes("glaberrima", "japonica", "PBM")
  ),
  filename = NULL, 
  fill = Set1[c(1:4)],
  main = expression("DE genes "*italic("vs")*". "*italic("japonica")*", PBM"),
  main.fontfamily = "Sans",
  cat.fontfamily = "Sans",
  fontfamily = "Sans",
  cat.fontface = "italic"
)
grid::grid.newpage()
grid::grid.draw(vd)
