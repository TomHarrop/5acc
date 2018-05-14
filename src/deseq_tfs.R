#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)

#############
# FUNCTIONS #
#############

AnnotateDESeqResultsTable <- function(results_set) {
  my_results <- data.table(as.data.frame(results_set), keep.rownames = TRUE)
  setnames(my_results, "rn", "MsuID")
  my_annotation <- my_results[, oryzr::LocToGeneName(unique(MsuID))]
  merge(my_results, my_annotation, by = "MsuID", all.x = TRUE)
}


PrintAnnotatedDESeqResultsTable <- function(annotated_results_table,
                                            file_name,
                                            outdir) {
  my_outfile <- file.path(outdir,
                          paste0(file_name, ".tab"))
  fwrite(annotated_results_table,
         my_outfile,
         sep = "\t")  
}  

###########
# GLOBALS #
###########

tfdb_file <- snakemake@input[["tfdb"]]
det_genes_file <- snakemake@input[["detected_genes"]]
dds_file <- snakemake@input[["dds"]]

alpha <- snakemake@params[["alpha"]]
lfcThreshold <- snakemake@params[["lfc_threshold"]]

cpus <- snakemake@threads[[1]]

log_file <- snakemake@log[["log"]]

all_outdir <- snakemake@params[["all_outdir"]]
sig_outdir <- snakemake@params[["sig_outdir"]]

# dev
# tfdb_file <- "output/010_data/tfdb.Rds"
# det_genes_file <- "output/060_tpm/detected_genes.Rds"
# dds_file <- "output/050_deseq/dds.Rds"
# cpus <- 8
# alpha <- 0.1
# lfcThreshold <- log(1.5, 2)
# outdir <- "output/deseq2/wald_tests/tf_only"

########
# MAIN #
########

BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read the dds
dds <- readRDS(dds_file)

# fix factors, thanks a lot DESeq
dds$continent <- factor(dds$continent)
dds$domestication <- factor(dds$domestication)

# get the list of detected genes
det_genes <- readRDS(det_genes_file)

# get a list of transcription factors
tfdb <- readRDS(tfdb_file)
all_tfs <- tfdb[, unique(`Protein ID`)]

# subset the dds for expressed TFs
kept_genes <- intersect(names(dds), intersect(all_tfs, det_genes))
dds_exp_tfs <- dds[kept_genes]

# remove japonica
dds_no_jp <- dds_exp_tfs[, dds_exp_tfs$accession != "japonica"]
dds_no_jp$accession <- droplevels(dds_no_jp$accession)

# DE analysis
dds_basal <- copy(dds_exp_tfs)
dds_stage_continent <- copy(dds_no_jp)
dds_wild_vs_dom_by_group <- copy(dds_no_jp)
dds_domestication <- copy(dds_no_jp)

design(dds_basal) <- ~ accession + stage
dds_basal <- DESeq(dds_basal, parallel = TRUE)

design(dds_stage_continent) <- ~ stage + continent + stage:continent
dds_stage_continent <- DESeq(dds_stage_continent, parallel = TRUE)

design(dds_wild_vs_dom_by_group) <- ~ stage + accession + stage:accession
dds_wild_vs_dom_by_group$accession <- 
  relevel(dds_wild_vs_dom_by_group$accession, "barthii") # barthii must be the base level
dds_wild_vs_dom_by_group <- DESeq(dds_wild_vs_dom_by_group, parallel = TRUE)

design(dds_domestication) <- ~ stage + domestication + stage:domestication
dds_domestication <- DESeq(dds_domestication, parallel = TRUE)

# basal branching TFs
results_basal <- results(dds_basal,
                         contrast = c("stage", "SM", "PBM"),
                         lfcThreshold = lfcThreshold,
                         alpha = alpha)

# de between stages is different for asian and african rice
results_stage_continent <- results(dds_stage_continent,
                                   name = "stageSM.continentAsia",
                                   alpha = alpha)

# accession-specific domestication differences
results_glab_barthii <- results(
  dds_wild_vs_dom_by_group,
  name = "stageSM.accessionglaberrima",
  alpha = alpha)

results_indica_rufipogon <- results(
  dds_wild_vs_dom_by_group,
  contrast = list("stageSM.accessionindica",
                  "stageSM.accessionrufipogon"),
  alpha = alpha)

# domestication / stage interaction
results_domestication <- results(
  dds_domestication,
  name = "stageSM.domesticationdomesticated",
  alpha = alpha)

# collect the results
sig_results <- c(
  "stage" = subset(results_basal,
                   padj < alpha),
  "interaction_stage_continent" = subset(results_stage_continent,
                                         padj < alpha),
  "interaction_stage_accession_indica" = subset(results_indica_rufipogon,
                                                padj < alpha),
  "interaction_stage_accession_glaberrima" = subset(results_glab_barthii,
                                                    padj < alpha),
  "domestication" = subset(results_domestication,
                           padj < alpha))

all_results <- c(
  "stage" = results_basal,
  "interaction_stage_continent" = results_stage_continent,
  "interaction_stage_accession_indica" = results_indica_rufipogon,
  "interaction_stage_accession_glaberrima" = results_glab_barthii,
  "domestication" = results_domestication)


# annotate
sig_results_annotated <- lapply(sig_results, AnnotateDESeqResultsTable)
all_results_annotated <- lapply(all_results, AnnotateDESeqResultsTable)

# write
sapply(names(sig_results_annotated), function(x)
  PrintAnnotatedDESeqResultsTable(sig_results_annotated[[x]],
                                  file_name = x,
                                  outdir = sig_outdir))
sapply(names(all_results_annotated), function(x)
  PrintAnnotatedDESeqResultsTable(all_results_annotated[[x]],
                                  file_name = x,
                                  outdir = all_outdir))

# write log
sessionInfo()
