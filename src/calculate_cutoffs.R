
library(data.table)
library(DESeq2)

###########
# GLOBALS #
###########


norm_counts_file <- snakemake@input[["norm_counts"]]
star_log_file <- snakemake@input[["star_logs"]]
bg_dds_file <- snakemake@input[["bg_dds"]]
feature_lengths_file <- snakemake@input[["feature_lengths"]]
tpm_file <- snakemake@input[["tpm"]]

log_file <- snakemake@log[["log"]]

detected_genes_file <- snakemake@output[["detected_genes"]]
tpm_with_calls_file <- snakemake@output[["tpm_with_calls"]]

# dev
# norm_counts_file <- "output/050_deseq/norm_counts.Rds"
# star_log_file <- "output/030_mapping/stats/star_logs.Rds"
# bg_dds_file <- "output/040_background-counts/dds_background.Rds"
# feature_lengths_file <- "output/010_data/feature_lengths.Rds"
# tpm_file <- "output/060_tpm/tpm.Rds"

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read data
star_log <- readRDS(star_log_file)
feature_lengths <- readRDS(feature_lengths_file)
norm_counts_wide <- readRDS(norm_counts_file)
dds_bg <- readRDS(bg_dds_file)
real_tpm <- readRDS(tpm_file)

# extract mu from STAR log
mu_table <- star_log[, .(library, mu = `Average mapped length`)]

# melt norm_counts
norm_counts <- melt(data.table(norm_counts_wide, keep.rownames = TRUE),
                    id.vars = "rn",
                    variable.name = "library",
                    value.name = "norm_counts")
setnames(norm_counts, "rn", "gene_name")


# normalize background counts
norm_counts_bg_wide <- counts(dds_bg, normalized = TRUE)
norm_counts_bg <- melt(data.table(norm_counts_bg_wide, keep.rownames = TRUE),
                       id.vars = "rn",
                       variable.name = "library",
                       value.name = "norm_counts")
setnames(norm_counts_bg, "rn", "gene_name")

# add feature lengths
norm_counts_fl <- merge(norm_counts, feature_lengths, by = "gene_name")
norm_counts_bg_fl <- merge(norm_counts_bg, feature_lengths, by = "gene_name")

# rename background genes
norm_counts_bg_fl[, gene_name := paste(gene_name,
                                       "background",
                                       sep = "_"),
                  by = gene_name]

# merge background and genic counts
combined_counts <- rbindlist(
  list(genic = norm_counts_fl, bg = norm_counts_bg_fl),
  idcol = "type")

# add mu
tpm_data <- merge(combined_counts, mu_table, by = "library")

# rename for tpm formula
setnames(tpm_data, c("mu", "norm_counts"), c("rl", "r.g"))

# calculate tpm
tpm_data[, T.g := r.g * rl / length, by = .(gene_name, library)]
tpm_data[, T.sum := sum(T.g), by = library]
tpm_data[, tpm := (r.g * rl * 1e6) / (length * T.sum)]

# calculate quantiles
cutoffs <- tpm_data[type == "bg",
         .(q95 = quantile(tpm, 0.95)),
         by = library]

# call genes
called_genes <- merge(tpm_data[type == "genic", .(library, gene_name, tpm)],
      cutoffs, by = "library")
 
called_genes[, lib_call := tpm > q95, by = .(library, gene_name)]

# merge calls only with real library TPM
called_tpm <- merge(real_tpm,
      called_genes[, .(library, gene_name, lib_call)],
      by = c("library", "gene_name"))

# calculate the real library cutoff
real_tpm_cutoffs <- called_tpm[lib_call == TRUE,
                               .(min_tpm = min(tpm)),
                               by = library]

# call at stage level 
called_tpm[, c("accession", "stage", "replicate") := tstrsplit(library, "_")]
called_tpm[, stage_call := sum(lib_call) >= 2,
           by = .(accession, stage, gene_name)]

# list detected genes
detected_genes <- called_tpm[stage_call == TRUE, unique(gene_name)]

# merge output
all_tpm_output <- merge(called_tpm, real_tpm_cutoffs, by = "library")
tpm_wide <- dcast(all_tpm_output, gene_name ~ library, value.var = "tpm")

# write output
saveRDS(all_tpm_output, tpm_with_calls_file)
saveRDS(detected_genes, detected_genes_file)

# write log
sessionInfo()
