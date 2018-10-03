#!/usr/bin/env Rscript

library(data.table)

star_logs <- fread("output/030_mapping/stats/star_logs.csv")

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

stage_order <- c("PBM" = "BM",
                 "SM" = "SM")

mapping_stats <- star_logs[, cols_to_keep, with = FALSE]
mapping_stats[, c("Species", "Stage", "Replicate") := tstrsplit(library, "_")]

mapping_stats[, Species := factor(plyr::revalue(Species, species_order),
                                  levels = species_order)]
mapping_stats[, Stage := factor(plyr::revalue(Stage, stage_order),
                                  levels = stage_order)]

setnames(mapping_stats, cols_to_keep, names(cols_to_keep))
mapping_stats[, library := NULL]

setcolorder(mapping_stats,
            c("Species", "Stage", "Replicate", names(cols_to_keep)[-1]))
setorder(mapping_stats, Species, Stage, Replicate)

fwrite(mapping_stats, "test/Table S1.csv")

