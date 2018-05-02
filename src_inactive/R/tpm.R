#!/usr/bin/Rscript

library(data.table)

################
# NOMENCLATURE #
################

# This script follows the naming in Wagner et al (2012; doi 
# 10.1007/s12064-012-0162-3). I don't know if this takes effective length of the
# region into account but the length of the mapped read is included in T.g.
# r.g = number of reads mapped to a particular gene region (counts)
#  fl = feature length
#  rl = read length (the average number of nucleotides mapped per read)
#   T = total number of transcripts sampled in a sequencing run

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

GenerateMessage("Calculate TPM")

# load feature.lengths
GenerateMessage("Loading feature lengths")
f.length.file <- "output/tpm/feature_lengths.Rds"
if (!file.exists(f.length.file)) {
  stop("feature_lengths.Rds not found")
}

feature.lengths <- data.table(readRDS(f.length.file), keep.rownames = TRUE)
setnames(feature.lengths, "rn", "gene")
setkey(feature.lengths, "gene")

# retrieve mu (average mapped length) from parsed STAR results
GenerateMessage("Reading average mapped length from STAR results")
star.log.file <- "output/mappingStats/starLogs.Rds"
if (!file.exists(star.log.file)) {
  stop("starLogs.Rds file not found")
}

star.logs <- readRDS(star.log.file)

# load DESeq2 results
GenerateMessage("Loading DESeq2 results")
dds.file <- "output/deseq2/dds.Rds"
if (!file.exists(dds.file)) {
  stop("dds.Rds file not found")
}

dds <- readRDS(dds.file)

# construct a long data.table of norm counts
normalized.counts.wide <- data.table(DESeq2::counts(dds, normalized = TRUE),
                                     keep.rownames = TRUE)
setnames(normalized.counts.wide, "rn", "gene")
normalized.counts <- reshape2::melt(
  normalized.counts.wide, id.vars = "gene", variable.name = "sample",
  value.name = "r.g")
setkey(normalized.counts, "gene")

# merge feature.lengths
GenerateMessage("Merging feature lengths")
tpm <- feature.lengths[normalized.counts]

# add mu per gene per sample
GenerateMessage("Merging mu (average mapped length)")
tpm[, rl := star.logs[Library == sample, `Average mapped length`], by = sample]
fl.nc <- copy(tpm)

# calculate T per gene
GenerateMessage("Calculating T")
tpm[, T.g := r.g * rl / Length, by = .(gene, sample)]

# sum T.g per sample
tpm[, T.sum := sum(T.g), by = sample]

# calculate TPM
GenerateMessage("Calculating TPM")
tpm[, tpm := (r.g * rl * 1e6) / (Length * T.sum)]

tpm.wide <- data.table(reshape2::dcast(tpm, gene ~ sample, value.var = "tpm"),
                       key = "gene")

# output results
out.dir <- "output/tpm"
GenerateMessage(paste("Saving output to", out.dir))
saveRDS(tpm[, .(gene, sample, tpm)], paste0(out.dir, "/tpm.Rds"))
saveRDS(tpm.wide, paste0(out.dir, "/tpm_wide.Rds"))
saveRDS(fl.nc, paste0(out.dir, "/fl_nc.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.tpm.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)
