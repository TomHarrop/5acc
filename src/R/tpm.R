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
GenerateMessageText <- function(message.text){
  paste0("[ ", date(), " ]: ", message.text)
}

message(GenerateMessageText("Calculate TPM"))

# load feature.lengths
message(GenerateMessageText("Loading feature lengths"))
f.length.file <- "output/tpm/feature_lengths.Rds"
if (!file.exists(f.length.file)) {
  stop("feature_lengths.Rds not found")
}

feature.lengths <- data.table(readRDS(f.length.file), keep.rownames = TRUE)
setnames(feature.lengths, "rn", "gene")
setkey(feature.lengths, "gene")

# retrieve mu (average mapped length) from parsed STAR results
message(GenerateMessageText("Reading average mapped length from STAR results"))
star.log.file <- "output/mappingStats/starLogs.Rds"
if (!file.exists(star.log.file)) {
  stop("starLogs.Rds file not found")
}

star.logs <- readRDS(star.log.file)

# load DESeq2 results
message(GenerateMessageText("Loading DESeq2 results"))
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
tpm <- feature.lengths[normalized.counts]

# add mu per gene per sample
tpm[, rl := star.logs[Library == sample, `Average mapped length`], by = sample]

# calculate T per gene
tpm[, T.g := r.g * rl / Length, by = .(gene, sample)]

# sum T.g per sample
tpm[, T.sum := sum(T.g), by = sample]

# calculate TPM
tpm[, tpm := (r.g * rl * 1e6) / (Length * T.sum)]


reshape2::dcast(tpm, gene ~ sample, value.var = "tpm")[1000:1005,]


############################################## 
# below is taken from Harold Pimentel's code #
##############################################
# calculate effective length per gene in each sample
# effective length = feature length - mu + 1
tpm[, li := Length - mu + 1, by = .(gene, sample)]

# calculate rate (counts per base)
tpm[, rate := Xi/li, by = .(gene, sample)]

# set the rate to zero if li == 0
tpm[li == 0, rate := 0]

# calculate the sum of rates per sample (Xj/lj)
tpm[, sum.of.rates := sum(rate), by = sample]

# calculate the TPM
tpm[, tpm := rate * 1e6 / sum.of.rates]


