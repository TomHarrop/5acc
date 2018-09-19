#!/usr/bin/env Rscript

library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}
GenerateMessage("Download ENSEMBL plant Compara and extract genewise homology")

# download ENSEMBL Plants Multi-species Compara file
GenerateMessage("Downloading file")
tmp <- tempfile(fileext = "tsv.gz")
compara.url <- 
  "ftp://ftp.ensemblgenomes.org/pub/plants/release-32/tsv/ensembl-compara/Compara.homologies.32.tsv.gz"
download.file(compara.url, tmp)

# read into R
GenerateMessage("Reading lines with zgrep and fread")
homology.table <- fread(paste("zgrep oryza", tmp))

# get column names
GenerateMessage("Reading column names")
names.table <- fread(paste("zcat", tmp, "| head -n 1"), verbose = FALSE)
setnames(homology.table, names(names.table))

# get homolog information for o. sativa
GenerateMessage("Extracting homolog information")
homology.species <- c("oryza_indica", "oryza_glaberrima", "oryza_barthii",
                      "oryza_rufipogon")
oryza.homology.table <- homology.table[
  species == "oryza_sativa" & homology_species %in% homology.species]

# get expected maximum homology identity by japonica gene
expected.identities <- oryza.homology.table[, .(
  expected.max.identity = max(identity),
  max.identity.homolog = homology_gene_stable_id[which.max(identity)]
), by = .(gene_stable_id, homology_species)]

# get RAP-MSU conversion table from oryzabase
GenerateMessage("Downloading RAP-MSU conversion table from oryzabase")
tmp <- tempfile(fileext = ".gz")
rapmsu.url <- "http://rapdb.dna.affrc.go.jp/download/archive/RAP-MSU.txt.gz"
download.file(rapmsu.url, tmp)
rapmsu.raw <- fread(paste("zcat", tmp), sep = "\t", header = FALSE,
                    col.names = c("gene_stable_id", "msu.tx.dirty"),
                    na.strings = "None")

# tidy up gene IDs
GenerateMessage("Tidying RAP-MSU conversion table")

rapmsu.raw[, gene_stable_id := toupper(gene_stable_id)]
TidyMsuId <- function(tx.dirty){
  tx.split <- unlist(strsplit(tx.dirty, ",", fixed = TRUE))
  genes <- gsub("\\.[[:digit:]]+$", "", tx.split)
  unique(genes[!is.na(genes)])
}
rapmsu <- rapmsu.raw[!is.na(gene_stable_id),
                     .(msu.id = TidyMsuId(msu.tx.dirty)),
                     by = gene_stable_id]
setkey(rapmsu, gene_stable_id)

# add msu IDs to expected.identities
GenerateMessage("Adding MSU IDs to homolog data")
setkey(expected.identities, gene_stable_id)
msu.homology <- merge(expected.identities, rapmsu)

# save output
out.dir <- "data/genome/ensembl"
GenerateMessage(paste("Saving output to", out.dir))
if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(homology.table, paste0(out.dir, "/Compara_homologies_oryza_raw.Rds"))
saveRDS(msu.homology, paste0(out.dir, "/msu_homology.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.genewise_homology_ensembl.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")

quit(save = "no", status = 0)
