#!/usr/bin/Rscript

library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

GenerateMessage("Calculate expression cutoffs")

# load real tpm
tpm.file <- "output/tpm/tpm.Rds"
GenerateMessage("Loading calculated TPM for genes")
if (!file.exists(tpm.file)) {
  stop("tpm.Rds not found")
}
tpm <- readRDS(tpm.file)

# load feature lengths and normalised counts for genic data
fl.nc.file <- "output/tpm/fl_nc.Rds"
GenerateMessage("Loading feature lengths and normalised counts for genes")
if (!file.exists(fl.nc.file)) {
  stop("fl_nc.Rds not found")
}
fl.nc <- readRDS(fl.nc.file)

# parse htseq-count output for intergenic regions
GenerateMessage("Reading count data for intergenic regions")
htseq.files <- list.files("output/shuffle/htseq", pattern = "^.*htseq-count$",
                          full.names = TRUE, recursive = TRUE)
names(htseq.files) <- gsub("^([[:alnum:]]+).*", "\\1", basename(htseq.files))
htseq.results <- lapply(htseq.files, read.table, header = FALSE, sep = "\t",
                        stringsAsFactors = FALSE, row.names = 1)
raw.intergenic.counts <- do.call(cbind, htseq.results)
colnames(raw.intergenic.counts) <- names(htseq.results)

# use DESeq2 to normalize intergenic counts
GenerateMessage("Normalizing intergenic counts with DESeq2")
meta.data <- data.frame(row.names = colnames(raw.intergenic.counts),
                        sample = colnames(raw.intergenic.counts))
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw.intergenic.counts[
    !grepl("^__", rownames(raw.intergenic.counts)), ],
  colData = meta.data,
  design = ~ 1
)
dds <- DESeq2::estimateSizeFactors(dds)
intergenic.counts.wide <- DESeq2::counts(dds, normalized = TRUE)
intergenic.counts <- data.table(reshape2::melt(intergenic.counts.wide,
                                               value.name =  "r.g"))
setnames(intergenic.counts, c("Var1", "Var2"), c("intergene", "sample"))
setkey(intergenic.counts, "intergene")

# merge gene counts and intergenic counts
GenerateMessage("Merging gene and intergenic counts")
feature.lengths <- unique(fl.nc[, .(gene, Length, rl)])
setkey(feature.lengths, "gene")
intergenic.counts <- feature.lengths[intergenic.counts]
intergenic.counts[, gene := paste("intergenic", as.character(gene), sep = ".")]
combined.counts <- rbind(fl.nc, intergenic.counts)

# calculate tpm on combined counts
GenerateMessage("Calculating tpm for genic and intergenic counts")
combined.counts[, T.g := r.g * rl / Length, by = .(gene, sample)]
combined.counts[, T.sum := sum(T.g), by = sample]
combined.counts[, tpm := (r.g * rl * 1e6) / (Length * T.sum)]
combined.counts[grep("^intergenic", gene), type := "intergenic"]
combined.counts[is.na(type), type := "genic"]

# calculate quantiles
GenerateMessage("Calculating threshhold for intergenic reads")
threshhold <- combined.counts[type == "intergenic",
                              .(q95 = quantile(tpm, 0.95)),
                              by = sample]
setkey(threshhold, "sample")

# call genes
GenerateMessage("Calling expressed genes")
genic.counts <- combined.counts[type == "genic", .(gene, sample, tpm)]
setkey(genic.counts, "sample")
genic.counts <- threshhold[genic.counts]
called.genes <- genic.counts[, .(call = tpm > q95), by = .(gene, sample)]
setkey(called.genes, "gene", "sample")

# merge real tpm
GenerateMessage("Merging calculated TPM with expression calls")
setkey(tpm, "gene", "sample")
tpm.with.calls <- called.genes[tpm]
tpm.with.calls[, species := substr(sample, 1, 1)]
tpm.with.calls[as.numeric(substr(sample, 2, 2)) <= 3,
               stage := "PBM"]
tpm.with.calls[as.numeric(substr(sample, 2, 2)) >= 4,
               stage := "SM"]
expressed.genes <- copy(tpm.with.calls)
detected.genes <- expressed.genes[, .(expressed.in.stage = sum(call) > 1),
                                  by = .(gene, species, stage)
                                  ][, .(exp = any(expressed.in.stage)),
                                    by = gene
                                    ][exp == TRUE, unique(gene)]

tpm.cutoffs <- tpm.with.calls[call == TRUE, .(min.tpm = min(tpm)), by = sample]

# save output
out.dir <- "output/tpm"
GenerateMessage(paste("Saving output to", out.dir))

saveRDS(tpm.with.calls, paste0(out.dir, "/tpm_with_calls.Rds"))
saveRDS(tpm.cutoffs, paste0(out.dir, "/tpm_cutoffs.Rds"))
saveRDS(detected.genes, paste0(out.dir, "/detected_genes.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.cutoffs.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)
