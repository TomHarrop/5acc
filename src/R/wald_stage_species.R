#!/usr/bin/env Rscript

library(DESeq2)
library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

GenerateMessage("Differential expression between stages within species")

# set up parallel processing
cpus <- as.numeric(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
if (is.na(cpus)) {
  cpus <- 1
}
GenerateMessage(paste("Allocating", cpus, "cpu(s)"))
BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# load data
GenerateMessage("Loading DESeq2 object")
dds.species.file <- "output/deseq2/wald_species/dds.species.Rds"
if (!file.exists(dds.species.file)) {
  stop("Couldn't find dds.species.Rds")
}
dds.species <- readRDS(dds.species.file)

# set up group variable and design
SummarizedExperiment::colData(dds.species)$group <- factor(paste(
  SummarizedExperiment::colData(dds.species)$stage,
  SummarizedExperiment::colData(dds.species)$accession, sep = "."))
DESeq2::design(dds.species) <- ~ group

# re-run DESeq2
GenerateMessage("Running filtered DE analysis with new design")
dds.stage.species <- DESeq2::DESeq(dds.species, parallel = TRUE)

# extract results
GenerateMessage("Extracting results")
ExtractResults <- function(accession, dds = dds.stage.species,
                           lfcThreshold = log(1.5, 2), alpha = 0.05) {
  message(paste0("                              ---> ", accession))
    contrast <- c(
    "group",
    paste("SM", accession, sep = "."),
    paste("PBM", accession, sep = "."))
  res <- DESeq2::results(
    object = dds, contrast = contrast,  lfcThreshold = lfcThreshold,
    alpha = alpha, parallel = TRUE)
  data.table(data.frame(res), keep.rownames = TRUE)
}
results.stage.species <- lapply(
  levels(SummarizedExperiment::colData(dds.stage.species)$accession),
  ExtractResults)
names(results.stage.species) <- levels(
  SummarizedExperiment::colData(dds.stage.species)$accession)
results.table <- rbindlist(results.stage.species, idcol = "accession")
setnames(results.table, "rn", "gene")
setkey(results.table, "accession", "gene")

# save output
GenerateMessage("Saving output")
out.dir <- "output/deseq2/wald_tests"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(dds.stage.species, paste0(out.dir, "/dds_stage_species.Rds"))
saveRDS(results.table, paste0(out.dir, "/stage_species_results_table.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.wald_stage_species.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)
