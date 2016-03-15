#!/usr/bin/env Rscript

library(DESeq2)
library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

GenerateMessage("Differential expression between species within stages")

# set up parallel processing
cpus <- as.numeric(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
if (is.na(cpus)) {
  cpus <- 1
}
GenerateMessage(paste("Allocating", cpus, "cpu(s)"))
BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# load data
GenerateMessage("Loading base DESeq2 object")
dds.file <- "output/deseq2/dds.Rds"
if (!file.exists(dds.file)) {
  stop("Couldn't find dds.Rds")
}
dds <- readRDS(dds.file)

# pre-filter based on tpm cutoffs
GenerateMessage("Removing undetected genes")
detected.genes.file <- "output/tpm/detected_genes.Rds"
if (!file.exists(detected.genes.file)) {
  stop("Couldn't find detected_genes.Rds")
}
detected.genes <- readRDS(detected.genes.file)

# re-run DESeq2 with expressed genes and new design
GenerateMessage("Running filtered DE analysis with new design")
design <- ~ stage + accession
dds.species <- DESeq2::DESeqDataSetFromMatrix(
  countData = DESeq2::counts(dds)[detected.genes, ],
  colData = SummarizedExperiment::colData(dds),
  design = design)
dds.species <- DESeq2::DESeq(dds.species, parallel = TRUE)

# extract results for each contrast
GenerateMessage("Extracting results")
accession <- levels(SummarizedExperiment::colData(dds.species)$accession)
combos <- combn(accession, m = 2, paste, collapse = ".")
contrasts <- lapply(combos, function(x) unlist(strsplit(x, split = "\\.")))
names(contrasts) <- combos

GetContrastResults <- function(x, dds = dds.species) {
  contrast <- c("accession", x[2], x[1])
  GenerateMessage(paste("--> Running contrast:", paste(contrast, collapse = ".")))
  DESeq2::results(
    dds.species,
    contrast = contrast,
    lfcThreshold = 1,
    alpha = 0.05,
    parallel = TRUE)
}
contrast.results.list <- lapply(
  contrasts, GetContrastResults, dds = dds.species)

# convert DESeq results to data.tables
GenerateMessage("Converting results to data.table and merging")
ConvertResultsToDt <- function(x) {
  x <- data.table(data.frame(x), keep.rownames = TRUE)
  setnames(x, "rn", "gene")
  x
}
results.tables <- lapply(contrast.results.list, ConvertResultsToDt)
contrast.results <- rbindlist(results.tables, idcol = TRUE)
setnames(contrast.results, ".id", "contrast")

# save output
GenerateMessage("Saving output")
out.dir <- "output/deseq2/wald_species"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(dds.species, paste0(out.dir, "/dds.species.Rds"))
saveRDS(contrast.results, paste0(out.dir, "/contrast_results.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.wald_species.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)



