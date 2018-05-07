#!/usr/bin/Rscript

library(data.table)

###########
# GLOBALS #
###########

feature_lengths_file <- snakemake@input[["feature_lengths"]]
star_log_file <- snakemake@input[["star_logs"]]
norm_counts_file <- snakemake@input[["norm_counts"]]


log_file <- snakemake@log[["log"]]

tpm_file <- snakemake@output[["tpm"]]
tpm_wide_file <- snakemake@output[["tpm_wide"]]
csv_file <- snakemake@output[["csv"]]

#########
# NOTES #
#########

# This script follows the naming in Wagner et al (2012; doi 
# 10.1007/s12064-012-0162-3). I don't know if this takes effective length of the
# region into account but the length of the mapped read is included in T.g.
# r.g = number of reads mapped to a particular gene region (counts)
#  fl = feature length
#  rl = read length (the average number of nucleotides mapped per read)
#   T = total number of transcripts sampled in a sequencing run

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# load data
norm_counts_df <- readRDS(norm_counts_file)
feature_lengths <- readRDS(feature_lengths_file)
star_log <- readRDS(star_log_file)

# extract mu from STAR log
mu_table <- star_log[, .(library, mu = `Average mapped length`)]

# get a long dt of normalised counts
norm_counts_wide <- data.table(norm_counts_df, keep.rownames = TRUE)
setnames(norm_counts_wide, "rn", "gene_name")
norm_counts <- melt(norm_counts_wide, id.vars = "gene_name",
                    variable.name = "library",
                    value.name = "r.g")

# add feature_lengths and mu
norm_counts_fl <- merge(norm_counts, feature_lengths, by = "gene_name")
tpm_data <- merge(norm_counts_fl, mu_table, by = "library")

# rename to match formula
setnames(tpm_data, c("mu"), c("rl"))

# calculate tpm
tpm_data[, T.g := r.g * rl / length, by = .(gene_name, library)]
tpm_data[, T.sum := sum(T.g), by = library]
tpm_data[, tpm := (r.g * rl * 1e6) / (length * T.sum)]

# extract data
tpm <- tpm_data[, .(library, gene_name, tpm)]
tpm_wide <- dcast(tpm, gene_name ~ library, value.var = "tpm")

# save output
saveRDS(tpm, tpm_file)
saveRDS(tpm_wide, tpm_wide_file)
fwrite(tpm, csv_file)

# write log
sessionInfo()
