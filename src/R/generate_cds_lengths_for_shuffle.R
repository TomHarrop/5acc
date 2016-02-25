#!/usr/bin/Rscript

library(data.table)

args <- commandArgs(TRUE)
outdir <- args[1]

# check for gtfLength
tpmDir <- "output/tpm"
if (!dir.exists(tpmDir)) {
  cat("tpmDir not found, exiting\n", file = stderr())
  quit(status = 1)
}

# import gtfLength
gtfLength <- data.table(readRDS(paste0(tpmDir, "/gtfLength.Rds")), keep.rownames = TRUE)

# import GTF
gtfFile <- "data/genome/os/Osativa_204_v7.0.gene_exons.cuffcomp.rRNAremoved.gtf"
gtf <- rtracklayer::import.gff(gtfFile, format = 'gtf', genome = 'Osativa_204_v7.0', feature.type="exon")

# make a dummy gff file with each 'gene' on the right chromosome but with the 
# coordinates 1 --> CDS length. This will be shuffled by bedtools shuffle so
# that the intergenic windows are the same size as acual genes.

# Uses by = rn to collect the duplicate seqids and strands that are unhelpfully
# returned from the GRanges object, then a second fast grouping to only keep the
# first of each.
dummyGff <- gtfLength[,.(
  seqid = unique(as.character(GenomeInfoDb::seqnames(gtf[gtf$gene_name == rn,]))),
  source = 'phytozomev10',
  type = 'CDS',
  start = 1,
  end = Length,
  score = ".",
  strand = unique(as.character(rtracklayer::strand(gtf[gtf$gene_name == rn,]))),
  phase = ".",
  attributes = paste0("ID=", rn)
), by = rn]

setkey(dummyGff, 'attributes')
dummyGff[, rn := NULL]

saveRDS(dummyGff, paste0(outdir, "/dummyGff.Rds"))

write.table(dummyGff, paste0(outdir, "/dummyGff.gff3"), sep = '\t', 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

