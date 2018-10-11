# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#!/usr/bin/env Rscript

library(rtracklayer)

###########
# GLOBALS #
###########

gff_file <- snakemake@input[["gff"]]
goi_list <- snakemake@params[["goi"]]
bed_file <- snakemake@output[["bed"]]

# dev
# gff_file <- "data/genome/os/Osativa_323_v7.0.gene_exons.gff3"
# goi_list <- list("LOC_Os01g04800", "LOC_Os10g33780")
# bed_file <- "test/goi.bed"

########
# MAIN #
########

# read the gff
os_exons <- import.gff3(gff_file, feature.type = "exon")

# search for genes in goi_list
goi <- unlist(goi_list)
ranges_list <- lapply(goi, function(x)
  os_exons[grep(x, os_exons$ID), ])

# combine results
all_ranges <- GRangesList(do.call(c, ranges_list))

# write output
export.bed(object = all_ranges, con = bed_file, format = "bed")

sessionInfo()