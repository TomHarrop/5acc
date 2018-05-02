#!/usr/bin/env Rscript

library(DESeq2)
library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

GenerateMessage("Differential expression between stages within continents")

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
  stop("Couldn't fine detected_genes.Rds")
}
detected.genes <- readRDS(detected.genes.file)

# set up group variable and design
col.data <- SummarizedExperiment::colData(dds)
col.data$group <- factor(paste(col.data$stage, col.data$continent, sep = "."))

# set up new DESeq2 object
dds.stage.continent <- DESeq2::DESeqDataSetFromMatrix(
  countData = DESeq2::counts(dds)[detected.genes, ],
  colData = col.data,
  design = ~ group)

# re-run DESeq2
GenerateMessage("Running filtered DE analysis with new design")
dds.stage.continent <- DESeq2::DESeq(dds.stage.continent, parallel = TRUE)

# extract results
GenerateMessage("Extracting results")
results.africa <- DESeq2::results(
  dds.stage.continent, contrast = c("group", "SM.Africa", "PBM.Africa"),
  lfcThreshold = log(1.5, 2), alpha = 0.05, parallel = TRUE)
results.asia <- DESeq2::results(
  dds.stage.continent, contrast = c("group", "SM.Asia", "PBM.Asia"),
  lfcThreshold = log(1.5, 2), alpha = 0.05, parallel = TRUE)

# convert to data.table
GenerateMessage("Converting results to data.table and merging")
results.africa.table <- data.table(data.frame(results.africa),
                                   keep.rownames = TRUE, key = "rn")
setnames(results.africa.table, "rn", "gene")
results.asia.table <- data.table(data.frame(results.asia),
                                   keep.rownames = TRUE, key = "rn")
setnames(results.asia.table, "rn", "gene")

# rbind tables
results.table <- rbindlist(
  list(asia = results.asia.table, africa = results.africa.table), idcol = TRUE)
setnames(results.table, ".id", "continent")

# save output
GenerateMessage("Saving output")
out.dir <- "output/deseq2/wald_tests"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(dds.stage.continent, paste0(out.dir, "/dds_stage_continent.Rds"))
saveRDS(results.table, paste0(out.dir, "/stage_continent_results_table.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.wald_stage_continent.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)


