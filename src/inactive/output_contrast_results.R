#!/usr/bin/env Rscript

library(data.table)

# function to annotate the tables
AnnotateResults <- function(x) {
  my.dt <- copy(x)
  my.annot <- my.dt[, oryzr::LocToGeneName(unique(gene))]
  merge(my.dt, my.annot, by.x = "gene", by.y = "MsuID")
}

# function to write output
WriteOutput <- function(x, file.name, out.dir = "explore/xlsx") {
  write.table(x, paste0(out.dir, "/", file.name), quote = FALSE, sep = "\t",
              na = "", row.names = FALSE)
}

# list of deseq results_tables
wald.files <- list.files("output/deseq2/wald_tests",
                         pattern = "results_table.Rds", full.names = TRUE,
                         include.dirs = FALSE)

# generate output names
names(wald.files) <- gsub("_results_table.Rds", "",
                          basename(wald.files), fixed = TRUE)

# read tables
wald.tables <- lapply(wald.files, readRDS)

# annotate results
wald.annotated <- lapply(wald.tables, AnnotateResults)

# extract sig genes
wald.annotated.sig <- lapply(wald.annotated, function(x) {
  x[padj < 0.05]
})

# save output
out.dir <- "output/deseq2/wald_tests/annotated_results"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
log.file <- paste0(out.dir, "/SessionInfo.txt")

lapply(names(wald.annotated), function(x)
  WriteOutput(wald.annotated[[x]],
              file.name = paste0(x, ".tab"),
              out.dir =  out.dir))

lapply(names(wald.annotated.sig), function(x)
  WriteOutput(wald.annotated.sig[[x]],
              file.name = paste0(x, "_sig.tab"),
              out.dir =  out.dir))

s.inf <- rutils::GitSessionInfo()
writeLines(s.inf, log.file)

quit(save = "no", status = 0)