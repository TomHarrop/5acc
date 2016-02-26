#!/usr/bin/Rscript

library(data.table)

args <- commandArgs(TRUE)
outdir <- args[1]

# messages
GenerateMessageText <- function(message.text){
  paste0("[ ", date(), " ]: ", message.text)
}

# check for gtfLength
tpmDir <- "output/tpm"
if (!dir.exists(tpmDir)) {
  cat("tpmDir not found, exiting\n", file = stderr())
  quit(status = 1)
}

# import gtfLength
message(GenerateMessageText("Loading feature.lengths"))
feature.lengths <- data.table(readRDS(paste0(tpmDir, "/feature_lengths.Rds")),
                              keep.rownames = TRUE, key = "rn")

# import GTF
gtf.file <-
  'data/genome/os/Osativa_323_v7.0.gene_exons.cuffcomp.rRNAremoved.gtf'
message(GenerateMessageText(paste("Reading features from file:",
                                  gtf.file, sep = "\n")))
gtf <- rtracklayer::import.gff(gtf.file, format = 'gtf',
                               genome = 'Osativa_323_v7.0', feature.type="exon")

# make a dummy gff file with each 'gene' on the right chromosome but with the 
# coordinates 1 --> CDS length. This will be shuffled by bedtools shuffle so
# that the intergenic windows are the same size as acual genes.

# make a data.table of gene_name, seqnames and strand and collapse
message(GenerateMessageText(
  "Extracting strand and chromosome information from gtf"))
gene.chromosome <- data.table(rn = gtf$gene_name,
                              seqid = as.character(GenomeInfoDb::seqnames(gtf)),
                              strand = as.character(rtracklayer::strand(gtf)),
                              key = "rn")
gene.chromosome <- unique(gene.chromosome)

# join gene.chromosome to feature.lengths and reformat to look like a gff
message(GenerateMessageText(
  "Joining chromosome information to feature lengths"))
dummy.gff <- gene.chromosome[feature.lengths, .(
  seqid,
  source = 'phytozomev10',
  type = 'CDS',
  start = 1,
  end = Length,
  score = ".",
  strand,
  phase = ".",
  attributes = paste0("ID=", rn)
)]
setkey(dummy.gff, 'attributes')

# save output
message(GenerateMessageText(paste("Saving dummyGff.gff3 to", outdir)))
saveRDS(dummy.gff, paste0(outdir, "/dummy_gff.Rds"))
write.table(dummy.gff, paste0(outdir, "/dummyGff.gff3"), sep = '\t', 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outdir, "/SessionInfo.generate_dummy_gff.txt")
writeLines(sInf, logLocation)

message(GenerateMessageText("Finished, exiting"))

quit(save = "no", status = 0)
