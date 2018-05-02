#!/usr/bin/env Rscript

library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}
GenerateMessage("Use biomaRt to get genes in GWAS regions")

# IRD network blocks the switch to GET after phytozome returns 302 from initial
# POST. Tell RCurl to use POST instead.
options(RCurlOptions = list(followlocation = TRUE, postredir = 2L))

# load crowell data
GenerateMessage("Loading GWAS regions")
gwas.crowell.binned.file <- "data/gwas/crowell/binned_gwas_crowell.Rds"
if (!file.exists(gwas.crowell.binned.file)) {
  stop("Couldn't find binned_gwas_crowell.Rds")
}
gwas.crowell.binned <- readRDS(gwas.crowell.binned.file)
#gwas.crowell.binned <- gwas.crowell.binned[1:10]

# setup phytozome
GenerateMessage("Setting up biomaRt")
phytozome <- biomaRt::useMart(biomart = "phytozome_mart",dataset = "phytozome",
                              host = "phytozome.jgi.doe.gov",
                              path = "/biomart/martservice/")

# setup query
bm.filters <- c("organism_id", "region_name", "region_start", "region_end")
bm.attributes <- c("gene_name1")

# query function
RetrieveRegionGenes <- function(region_name, region_start, region_end) {
  bm.values <- list(
    organism_id = "323", # O. sativa
    region_name = region_name,
    region_start = region_start,
    region_end = region_end)
  unique(biomaRt::getBM(bm.attributes, bm.filters, bm.values,
                        phytozome, uniqueRows = TRUE)$gene_name1)
}

# format query
gwas.crowell.binned[, c(
  "region_name", "region_start", "region_end") :=
    tstrsplit(Bin_id, split = "[^[:alnum:]]+"),
  by = Bin_id]
gwas.crowell.binned[, region_name := paste0("Chr", CHROMOSOME)]

# run query
n.bins <- gwas.crowell.binned[, length(unique(Bin_id))]
GenerateMessage(paste0("Querying phytozome for ", n.bins, " regions"))
bin.to.gene <- gwas.crowell.binned[, .(
  gene = RetrieveRegionGenes(region_name, region_start, region_end)),
  by = Bin_id]

# merge results
GenerateMessage("Merging results")
setkey(gwas.crowell.binned, Bin_id)
setkey(bin.to.gene, Bin_id)
gwas.crowell.genes <- bin.to.gene[gwas.crowell.binned, allow.cartesian = TRUE]

# save output
out.dir <- "output/gwas/crowell"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(gwas.crowell.genes, paste0(out.dir, "/gwas_crowell_genes.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.crowell_genes.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)
