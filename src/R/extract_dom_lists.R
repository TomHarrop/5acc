#!/usr/bin/env Rscript

library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

GenerateMessage(
  "Generate tables of wald test results for domestication contrasts")

# data files
GenerateMessage("Loading data")
dom.file <- "output/deseq2/wald_tests/domestication_results_table.Rds"
dom.asia.file <-
  "output/deseq2/wald_tests/domestication_asia_results_table.Rds"
dom.continent.file <- 
  "output/deseq2/wald_tests/domestication_by_continent_results_table.Rds"
stage.results.file <-
  "output/deseq2/wald_tests/stage_species_results_table.Rds"
for (dom.file in list(
  dom.all.file, dom.asia.file, dom.continent.file, stage.results.file)) {
  if (!file.exists(dom.file)) {
    stop(paste("Couldn't find", basename(dom.file)))
    quit(save = "no", status = "1")
  }
}
dom <- readRDS(dom.all.file)
dom.asia <- readRDS(dom.asia.file)
dom.continent <- readRDS(dom.continent.file)
stage.results.table <- readRDS(stage.results.file)

# join asia and continent results and do some munging
GenerateMessage("Joining and casting data.tables")
dom[, domestication := "all"]
dom.all.long <- rbind(dom.asia, dom.continent, dom)
dom.padj <- dcast(dom.all.long, gene ~ domestication, value.var = "padj")

# add stage results
GenerateMessage("Adding stage Wald test results")
setkey(dom.padj, gene)
setkey(stage.results.table, gene, accession)
stage.l2fc.dom.padj <- dom.padj[stage.results.table]
old.names <- c("all", "africa", "asia", "indica", "japonica")
setnames(stage.l2fc.dom.padj, old.names, paste("dom", old.names, sep = "_"))

# save output
out.dir <- "output/deseq2/wald_tests"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(stage.l2fc.dom.padj, paste0(out.dir, "/stage_l2fc_dom_padj.Rds"))
saveRDS(dom.all.long, paste0(out.dir, "/dom_all_long.Rds"))
saveRDS(dom.padj, paste0(out.dir, "/dom_padj.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.extract_dom_lists.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)
