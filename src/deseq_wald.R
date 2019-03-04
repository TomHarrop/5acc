#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)

#############
# FUNCTIONS #
#############

AnnotateDESeqResultsTable <- function(results_set, full_annotation) {
  merge(results_set,
        full_annotation,
        by = "gene_id",
        all.x = TRUE,
        all.y = FALSE)
}

DESeqResultsToDataTable <- function(res, dds) {
  my_res <- data.table(as.data.frame(res), keep.rownames = TRUE)
  setnames(my_res, "rn", "gene_id")
  my_res[, design := as.character(design(dds))[[2]]]
  my_res[, wald_test := sub("Wald test p-value: ", "", mcols(res)[5, 2])]
  return(my_res)
}

GetAccessionContrastResults <- function(contrast, name, dds) {
  my_contrast <- c(name, contrast[[2]], contrast[[1]])
  DESeqResultsToDataTable(results(dds,
                                  contrast = my_contrast,
                                  lfcThreshold = lfcThreshold,
                                  alpha = alpha),
                          dds)
}

PrintAnnotatedDESeqResultsTable <- function(annotated_results_table,
                                            file_name,
                                            all_outdir,
                                            sig_outdir) {
  my_allfile <- file.path(all_outdir,
                          paste0(file_name, ".csv"))
  my_sigfile <- file.path(sig_outdir,
                          paste0(file_name, ".csv"))
  fwrite(annotated_results_table,
         my_allfile)  
  fwrite(annotated_results_table[padj < alpha],
         my_sigfile)  
}  

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds"]]

alpha <- snakemake@params[["alpha"]]
lfcThreshold <- snakemake@params[["lfc_threshold"]]

all_outdir <- snakemake@params[["all_outdir"]]
sig_outdir <- snakemake@params[["sig_outdir"]]

cpus <- snakemake@threads[[1]]

log_file <- snakemake@log[["log"]]

# dev
# dds_file <- "output/050_deseq/filtered_dds.Rds"
# alpha <- 0.1
# lfcThreshold <- log(1.5, 2)

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

# remove japonica
dds_no_jp <- dds[, dds$accession != "japonica"]
dds_no_jp$accession <- droplevels(dds_no_jp$accession)

# DE TESTS #

# domestication : stage interaction (no japonica)
# +ve L2FC indicates increased L2FC between PBM and SM in domesticated accessions
dds_domestication <- copy(dds_no_jp)
design(dds_domestication) <- ~ domestication + stage + domestication:stage
dds_domestication <- DESeq(dds_domestication,
                           parallel = TRUE)
res_dds_domestication <- DESeqResultsToDataTable(
  results(dds_domestication,
          name = "domesticationdomesticated.stageSM",
          alpha = alpha),
  dds_domestication)

# compare interaction in indica, japonica and glaberrima separately,
# each compared to their wild relative
# covers wald_domestication_asia, wald_stage_accession, wald_domestication_continent
# formerly "domestication asia"
# rufipogon must be base level
# +ve L2FC indicates increased L2FC between PBM and SM in domesticated accession
dds_stage_accession <- copy(dds) 
dds_stage_accession$accession <- relevel(dds_stage_accession$accession,
                                         "rufipogon") 
design(dds_stage_accession) <- ~ stage + accession + stage:accession
dds_stage_accession <- DESeq(dds_stage_accession,
                             parallel = TRUE)
res_dds_stage_accession_japonica <- DESeqResultsToDataTable(
  results(dds_stage_accession,
          name = "stageSM.accessionjaponica",
          alpha = alpha),
  dds_stage_accession)
res_dds_stage_accession_indica <- DESeqResultsToDataTable(
  results(dds_stage_accession,
          name = "stageSM.accessionindica",
          alpha = alpha),
  dds_stage_accession)
res_dds_stage_accession_glaberrima <- DESeqResultsToDataTable(
  results(dds_stage_accession,
          contrast = list("stageSM.accessionglaberrima",
                          "stageSM.accessionbarthii"),
          alpha = alpha),
  dds_stage_accession)

# stage controlling for accession OR accession controlling for stage
dds_accession <- copy(dds)
design(dds_accession) <- ~ stage + accession # same model for dds_stage
dds_accession <- DESeq(dds_accession,
                       parallel = TRUE)
res_dds_accession_stage <- DESeqResultsToDataTable(
  results(dds_accession,
          name = "stage_SM_vs_PBM",
          alpha = alpha,
          lfcThreshold = lfcThreshold),
  dds_accession)

accessions <- levels(colData(dds_accession)$accession)
combos <- combn(accessions, m = 2, paste, collapse = "_")
contrasts <- lapply(combos, function(x) unlist(strsplit(x, split = "_")))
names(contrasts) <- combos
combo_results <- lapply(contrasts,
                        GetAccessionContrastResults,
                        name = "accession",
                        dds = dds_accession)
res_dds_accession <- rbindlist(combo_results)

# DE between stages within species OR DE between accessions within stage
dds_stage_between_species <- copy(dds)
dds_stage_between_species$group <- factor(
  paste(dds_stage_between_species$accession,
        dds_stage_between_species$stage,
        sep = "_"))
design(dds_stage_between_species) <- ~ group # same model for between stage within species
dds_stage_between_species <- DESeq(dds_stage_between_species,
                                   parallel = TRUE)
accessions <- levels(colData(dds_stage_between_species)$accession)
stages <- levels(colData(dds_stage_between_species)$stage)
acc_combos <- combn(accessions, m = 2, paste, collapse = "_")
acc_contrasts_1 <- lapply(acc_combos, function(x)
  paste(unlist(strsplit(x, split = "_")), stages[[1]], sep = "_"))
acc_contrasts_2 <- lapply(acc_combos, function(x)
  paste(unlist(strsplit(x, split = "_")), stages[[2]], sep = "_"))
stage_contrasts <- zipup(paste(accessions, stages[[1]], sep = "_"),
                         paste(accessions, stages[[2]], sep = "_"))
acc_combo_results <- lapply(c(acc_contrasts_1, acc_contrasts_2),
                            GetAccessionContrastResults,
                            name = "group",
                            dds = dds_stage_between_species)
stage_combo_results <- lapply(stage_contrasts,
                              GetAccessionContrastResults,
                              name = "group",
                              dds = dds_stage_between_species)
res_dds_stage_between_species <- rbindlist(acc_combo_results)
res_dds_stage_within_species <- rbindlist(stage_combo_results)

# DE between stages within asian OR african accessions (no japonica)
dds_stage_continent <- copy(dds_no_jp)
dds_stage_continent$group <- factor(
  paste(dds_stage_continent$stage,
        dds_stage_continent$continent,
        sep = "_"))
design(dds_stage_continent) <- ~ group
dds_stage_continent <- DESeq(dds_stage_continent,
                             parallel = TRUE)
res_dds_stage_continent_africa <- DESeqResultsToDataTable(
  results(dds_stage_continent,
          contrast = c("group",
                       "SM_Africa",
                       "PBM_Africa"),
          alpha = alpha),
  dds_stage_continent)
res_dds_stage_continent_asia <- DESeqResultsToDataTable(
  results(dds_stage_continent,
          contrast = c("group",
                       "SM_Asia",
                       "PBM_Asia"),
          alpha = alpha),
  dds_stage_continent)

# make a list of all the results
all_wald_results <- list(
  "domestication" = res_dds_domestication,
  "stage_accession_japonica" = res_dds_stage_accession_japonica,
  "stage_accession_indica" = res_dds_stage_accession_indica,
  "stage_accession_glaberrima" = res_dds_stage_accession_glaberrima,
  "stage_between_species" = res_dds_stage_between_species,
  "stage_within_species" = res_dds_stage_within_species,
  "stage_continent_africa" = res_dds_stage_continent_africa,
  "stage_continent_asia" = res_dds_stage_continent_asia,
  "accession" = res_dds_accession,
  "stage" = res_dds_accession_stage)

# generate an annotation
full_annotation <- oryzr::LocToGeneName(unique(rownames(dds)))
setnames(full_annotation, "MsuID", "gene_id")

# annotate the results
annotated_wald_results <- lapply(all_wald_results,
                                 AnnotateDESeqResultsTable,
                                 full_annotation = full_annotation)

# print the results
lapply(names(annotated_wald_results), function(x)
  PrintAnnotatedDESeqResultsTable(annotated_wald_results[[x]],
                                  x,
                                  all_outdir,
                                  sig_outdir))

# write log
sessionInfo()
