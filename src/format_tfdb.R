#!/usr/bin/env Rscript

library(data.table)

###########
# GLOBALS #
###########

raw_tfdb_file <- snakemake@input[["tfdb_file"]]

log_file <- snakemake@log[["log"]]

formatted_tfdb_file <- snakemake@output[["tfdb_file"]]

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read tfdb
tfdb_raw <- fread(raw_tfdb_file)

# remove duplicates
unique_tfs <- unique(tfdb_raw, by = c("Species", "Family", "Protein ID"))

# extract OS TFs
os <- unique(unique_tfs[Species == 'Oryza sativa subsp. japonica', .(
  "Protein ID" = gsub("\\..*", '', `Protein ID`),
  Family
)])

# write output
saveRDS(os, formatted_tfdb_file)

# write log
sessionInfo()

