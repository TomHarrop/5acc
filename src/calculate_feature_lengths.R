#!/usr/bin/env Rscript

library(data.table)
library(rtracklayer)

###########
# GLOBALS #
###########

os_gtf_file <- snakemake@input[["os_gtf_file"]]
log_file <- snakemake@log[["log"]]
feature_lengths_file <- snakemake@output[["feature_lengths"]]

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# import GTF
os_gff_exons <- import.gff(os_gtf_file,
                  format = 'gtf',
                  genome = 'Osativa_323_v7.0',
                  feature.type="exon")

# merge overlapping exons per-gene
grl <- GenomicRanges::reduce(split(os_gff_exons,
                                   elementMetadata(os_gff_exons)$gene_name))
gtf_reduced <- unlist(os_gff_exons, use.names = FALSE)

# add metadata
elementMetadata(gtf_reduced)$widths <- width(gtf_reduced)

# calculate feature lengths with dplyr
feature_lengths <- as.data.table(gtf_reduced)[
  , .(length = sum(width)), by = gene_name]

# write output
saveRDS(feature_lengths, feature_lengths_file)

# write log
sessionInfo()

