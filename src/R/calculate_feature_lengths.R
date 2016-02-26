#!/usr/bin/Rscript

library(rtracklayer)
library(dplyr)

# messages
GenerateMessageText <- function(message.text){
  paste0("[ ", date(), " ]: ", message.text)
}

message(GenerateMessageText("Calculate feature lengths from GFF3"))

# code from devon ryan:
# http://seqanswers.com/forums/showpost.php?p=129175&postcount=3

# import GTF
gtf.file <-
  'data/genome/os/Osativa_323_v7.0.gene_exons.cuffcomp.rRNAremoved.gtf'
message(GenerateMessageText(paste("Reading features from file:",
                                  gtf.file, sep = "\n")))
gtf <- import.gff(gtf.file, format = 'gtf', genome = 'Osativa_323_v7.0',
                  feature.type="exon")

# reduce ranges by gene_name (MSU ID), i.e. merge overlapping exons
message(GenerateMessageText("Merging exons"))
grl <- reduce(split(gtf, elementMetadata(gtf)$gene_name))
gtf.reduced <- unlist(grl, use.names = FALSE)

# add metadata
elementMetadata(gtf.reduced)$gene_name <- rep(names(grl), elementLengths(grl))
elementMetadata(gtf.reduced)$widths <- width(gtf.reduced)

# calculate feature lengths with dplyr
message(GenerateMessageText(
  "Grouping merged exons and calculating gene lengths"))
output <- group_by(as.data.frame(gtf.reduced), gene_name) %>%
  summarize(length = sum(widths))
feature.lengths <- data.frame(Length = output$length,
                              row.names = output$gene_name)

# save output
outDir <- "output/tpm"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

message(GenerateMessageText(paste("Saving feature_lengths.Rds to", outDir)))
if (file.exists(paste0(outDir, "/feature_lengths.Rds"))){
  warning("Replacing existing file")
}
saveRDS(feature.lengths, paste0(outDir, "/feature_lengths.Rds"))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(outDir, "/SessionInfo.calculate_feature_lengths.txt")
writeLines(sInf, logLocation)

message(GenerateMessageText("Finished, exiting"))

quit(save = "no", status = 0)
