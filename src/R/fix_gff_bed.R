#!/usr/bin/Rscript

args <- commandArgs(TRUE)
outdir <- args[1]

bed9to6 <- function(bed9file) {
  bed9 <- read.table(bed9file, sep = "\t", stringsAsFactors = FALSE)
  bed9[, 5] <- 0
  return(bed9[, 1:6])
}

bed9files <- list.files(outdir, pattern = "bed$", full.names = TRUE)

bed6 <- lapply(bed9files, bed9to6)
names(bed6) <- paste0(gsub("\\.bed.*", "",
                           sapply(bed9files, basename)), '.bed6')

bed6[["irgsp1_rRNA_tRNA.bed6"]]$V1 <- gsub(
  "0(\\d)", "\\1", gsub("chr", "Chr", bed6[["irgsp1_rRNA_tRNA.bed6"]]$V1,
                        fixed = TRUE))

for (i in 1:length(bed6)) {
  outpath <- paste0(outdir, "/", names(bed6)[i])
  write.table(bed6[i], file = outpath, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
}