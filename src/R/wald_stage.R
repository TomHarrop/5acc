#!/usr/bin/env Rscript

library(DESeq2)
library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

GenerateMessage("Differential expression between stages across species")

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

# extract results for stage
GenerateMessage("Extracting results")
stage.results <- DESeq2::results(
  dds.species, contrast = c("stage", "SM", "PBM"), lfcThreshold = log(1.5, 2),
  alpha = 0.05, parallel = TRUE)
stage.results.table <- data.table(data.frame(stage.results),
                                  keep.rownames = TRUE, key = "rn")
setnames(stage.results.table, "rn", "gene")

# save output
GenerateMessage("Saving output")
out.dir <- "output/deseq2/wald_tests"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(stage.results.table, paste0(out.dir, "/stage_results_table.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.wald_stage.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)