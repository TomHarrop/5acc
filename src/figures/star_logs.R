#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# Table S5: mapping stats

library(data.table)

star_logs <- fread(snakemake@input[["star_logs"]])

cols_to_keep <- c("library" = "library",
                  "Input reads (M)" = "Number of input reads",
                  "Uniquely mapped reads (M)" = "Uniquely mapped reads number",
                  "Uniquely mapped reads in genes (M)" = "Number of reads in genes",
                  "Unique mapping %" = "Uniquely mapped reads %")

species_order <- c("or" = "Oryza rufipogon",
                   "osi" = "Oryza sativa indica",
                   "osj" = "Oryza sativa japonica",
                   "ob" = "Oryza barthii",
                   "og" = "Oryza glaberrima")

stage_order <- c("PBM" = "IM",
                 "SM" = "DM")

# select columns
mapping_stats <- star_logs[, cols_to_keep, with = FALSE]
mapping_stats[, c("Species", "Stage", "Replicate") := tstrsplit(library, "_")]

# divide by 1 million
mapping_stats[, cols_to_keep[2:4] := lapply(.SD, function(x) round(x/1e6, 1)),
              .SDcols = cols_to_keep[2:4]]

# reorder
mapping_stats[, Species := factor(plyr::revalue(Species, species_order),
                                  levels = species_order)]
mapping_stats[, Stage := factor(plyr::revalue(Stage, stage_order),
                                  levels = stage_order)]

# rename
setnames(mapping_stats, cols_to_keep, names(cols_to_keep))
mapping_stats[, library := NULL]

# order cols
setcolorder(mapping_stats,
            c("Species", "Stage", "Replicate", names(cols_to_keep)[-1]))
setorder(mapping_stats, Species, Stage, Replicate)

# write output
fwrite(mapping_stats[Species != "Oryza sativa japonica"],
       snakemake@output[["table1"]])

sessionInfo()
