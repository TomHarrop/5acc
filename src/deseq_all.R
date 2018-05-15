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
cpus <- 8
dds_file <- "output/050_deseq/filtered_dds.Rds"
alpha <- 0.1
lfcThreshold <- log(1.5, 2)


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

# copy dds for different models
dds_domestication <- copy(dds_no_jp)
design(dds_domestication) <- ~ domestication + stage + domestication:stage

dds_stage_accession <- copy(dds_no_jp) # formerly "domestication asia"
design(dds_stage_accession) <- ~ stage + accession + stage:accession

dds_domestication_continent <- copy(dds_no_jp)
dds_domestication_continent$group <- factor(
  paste(dds_domestication_continent$continent,
        dds_domestication_continent$domestication,
        sep = "_"))
design(dds_domestication_continent) <- ~ group + stage + group:stage

dds_accession <- copy(dds_no_jp)
design(dds_accession) <- ~ stage + accession # same model for dds_stage

dds_stage_between_species <- copy(dds_no_jp)
dds_stage_between_species$group <- factor(
  paste(dds_stage_between_species$accession,
        dds_stage_between_species$stage,
        sep = "_"))
model(dds_stage_between_species) <- ~ group # same model for between stage within species

dds_stage_continent <- copy(dds_no_jp)
dds_stage_continent$group <- factor(
  paste(dds_stage_continent$stage,
        dds_stage_continent$continent,
        sep = "_"))
model(dds_stage_continent) <- ~ group




