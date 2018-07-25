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

# manually add ALOG genes
ALOG <- c('LOC_Os07g04670', 'LOC_Os02g07030', 'LOC_Os06g46030',
          'LOC_Os02g41460', 'LOC_Os04g43580', 'LOC_Os10g33780',
          'LOC_Os02g56610', 'LOC_Os01g61310', 'LOC_Os05g39500',
          'LOC_Os05g28040')
os_with_alog <- rbind(os, 
                      data.table(`Protein ID` = ALOG,
                                 Family = "ALOG"))

# write output
saveRDS(os_with_alog, formatted_tfdb_file)

# write log
sessionInfo()

