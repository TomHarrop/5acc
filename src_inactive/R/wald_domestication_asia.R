#!/usr/bin/env Rscript

library(DESeq2)
library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}
GenerateMessage(
  "Calculate domestication effect separately for japonica and indica")

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

# setup new deseq obect
col.data <- SummarizedExperiment::colData(dds)
col.data$accession <- factor(
  col.data$accession,
  levels = c("rufipogon", "indica", "japonica", "barthii", "glaberrima"))

dds.domestication.asia <- DESeq2::DESeqDataSetFromMatrix(
  countData = DESeq2::counts(dds)[detected.genes, ],
  colData = col.data,
  design = ~ stage + accession + stage:accession)

# rerun with updated design
GenerateMessage("Running DESeq2 with updated design")
dds.domestication.asia <- DESeq2::DESeq(dds.domestication.asia,
                                        parallel = TRUE)

# extract results
GenerateMessage("Extracting results")
results.japonica <- DESeq2::results(
  dds.domestication.asia,
  name = "stageSM.accessionjaponica",
  alpha = 0.05, parallel = TRUE)
results.indica <- DESeq2::results(
  dds.domestication.asia,
  name = "stageSM.accessionindica",
  alpha = 0.05, parallel = TRUE)

# combine results
GenerateMessage("Combining results")
results.list <- list(japonica = results.japonica, indica = results.indica)
ConvertResultsToDt <- function(x) {
  x <- data.table(data.frame(x), keep.rownames = TRUE, key = "rn")
  setnames(x, "rn", "gene")
}
results.tables <- lapply(results.list, ConvertResultsToDt)
results.table <- rbindlist(results.tables, idcol = TRUE)
setnames(results.table, ".id", "domestication")

# find genes that are significant twice: THIS WORKS
# results.table.wide <- dcast(
#   results.table, gene ~ domestication, value.var = c("log2FoldChange", "padj"))
# results.table.wide[padj_indica < 0.05 & padj_japonica < 0.05, unique(gene)]

# save output
GenerateMessage("Saving output")
out.dir <- "output/deseq2/wald_tests"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(dds.domestication.asia,
        paste0(out.dir, "/dds_domestication_asia.Rds"))
saveRDS(results.table, paste0(out.dir, "/domestication_asia_results_table.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir,
                      "/SessionInfo.wald_domestication_asia.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)


