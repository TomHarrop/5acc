#!/usr/bin/env Rscript

library(Rsamtools)

###########
# GLOBALS #
###########

rrna_bedfile <- "data/genome/os/rRna.combined.bed"
bamfile_dir <- "output/osj/STAR"
cpus <- 8

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# set up multiprocessing
BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# prepare annotations
features <- rtracklayer::import.bed(rrna_bedfile)

# find bamfiles
bamfile_paths <- list.files(bamfile_dir,
           recursive = FALSE,
           full.names = TRUE,
           pattern = "Aligned.out.bam")

bamfiles <- BamFileList(bamfile_paths)

# run counts
counts_ga <- GenomicAlignments::summarizeOverlaps(
  features = features,
  reads = bamfiles,
  mode = "Union",
  singleEnd = FALSE,
  ignore.strand = TRUE,
  fragments = TRUE)
