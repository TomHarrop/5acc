#!/usr/bin/env Rscript

library(DESeq2)
library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}
GenerateMessage("Calculate interaction between domestication and stage")

# set up parallel processing
cpus <- as.numeric(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
if (is.na(cpus)) {
  cpus <- 1
}
GenerateMessage(paste("Allocating", cpus, "cpu(s)"))
BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# new design
design <- ~ group + stage + group:stage

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
  stop("Couldn't fine detected_genes.Rds")
}
detected.genes <- readRDS(detected.genes.file)

# set up grouping variable
colData <- SummarizedExperiment::colData(dds)
colData$group <- factor(paste0(colData$continent, colData$domestication))

# set up new DESeq2 object
dds.domestication.by.continent <- DESeq2::DESeqDataSetFromMatrix(
  countData = DESeq2::counts(dds)[detected.genes, ],
  colData = colData,
  design = design)

# run DESeq2
GenerateMessage("Running DESeq2 with updated design")
dds.domestication.by.continent <- DESeq2::DESeq(dds.domestication.by.continent,
                                                parallel = TRUE)

# extract results
GenerateMessage("Extracting results")
results.asia.contrast <- list("groupAsiawild.stageSM",
                              "groupAsiadomesticated.stageSM")
results.asia <- DESeq2::results(
  dds.domestication.by.continent,
  contrast = results.asia.contrast,
  alpha = 0.05, parallel = TRUE)
results.africa <- DESeq2::results(
  dds.domestication.by.continent,
  name = "groupAfricawild.stageSM",
  alpha = 0.05, parallel = TRUE)

# combine results
GenerateMessage("Combining results")
results.list <- list(asia = results.asia, africa = results.africa)

ConvertResultsToDt <- function(x) {
  x <- data.table(data.frame(x), keep.rownames = TRUE, key = "rn")
  setnames(x, "rn", "gene")
}

results.tables <- lapply(results.list, ConvertResultsToDt)
results.table <- rbindlist(results.tables, idcol = TRUE)
setnames(results.table, ".id", "domestication")

# save output
GenerateMessage("Saving output")
out.dir <- "output/deseq2/wald_domestication_by_continent"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(dds.domestication.by.continent, paste0(out.dir, "/dds_domestication_by_continent.Rds"))
saveRDS(results.table, paste0(out.dir, "/results_table.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.wald_domestication_by_continent.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)

