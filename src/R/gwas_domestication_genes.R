#!/usr/bin/env Rscript

library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}
GenerateMessage("Look for domestication-associated genes in GWAS regions")

# load data
GenerateMessage("Loading data")
gwas.crowell.genes.file <- "output/gwas/crowell/gwas_crowell_genes.Rds"
stage.l2fc.dom.padj.file <- 
  "output/deseq2/wald_tests/stage_l2fc_dom_padj.Rds"

for (data.file in list(
  gwas.crowell.genes.file, stage.l2fc.dom.padj.file
)) {
  if (!file.exists(data.file)) {
    stop(paste("Couldn't find", basename(data.file)))
    quit(save = "no", status = "1")
  }
}
gwas.crowell.genes <- readRDS(gwas.crowell.genes.file)[!is.na(gene)]
stage.l2fc.dom.padj <- readRDS(stage.l2fc.dom.padj.file)

# join domestication to gwas
GenerateMessage("Joining domestication and gwas results")
setkey(stage.l2fc.dom.padj, gene, accession)
setkey(gwas.crowell.genes, gene)
gwas.with.p <- stage.l2fc.dom.padj[gwas.crowell.genes]

# collapse genes by TRAIT and SUB_POP
GenerateMessage("Collapsing trait and subpop to single rows")
gwas.with.p[, trait.collapsed := paste(unique(TRAIT), collapse = ", "),
            by = gene]
gwas.with.p[, subpop.collapsed := paste(unique(SUB_POP), collapse = ", "),
            by = gene]

# select rows
GenerateMessage("Tidying up the table")
setkey(gwas.with.p, gene, accession)
gwas.with.p <- unique(gwas.with.p)[, .(
  gene, dom_all, dom_africa, dom_asia, dom_japonica, dom_indica, accession,
  baseMean, log2FoldChange, lfcSE, padj, Bin_id, trait.collapsed,
  subpop.collapsed, region_name, region_start, region_end
)]

# save output
out.dir <- "output/gwas"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(gwas.with.p, paste0(out.dir, "/gwas_domestication.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.domestication_genes.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)

