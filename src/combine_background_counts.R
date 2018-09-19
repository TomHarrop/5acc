#!/usr/bin/Rscript

library(data.table)
library(DESeq2)

#############
# FUNCTIONS #
#############

GenerateName <- function(x, file) {
  paste(
    gsub("^.+/", "", dirname(x)),
    sub(file, "", basename(x)),
    sep = "_")
}

###########
# GLOBALS #
###########

count_file_list <- snakemake@input[["count_files"]]
counts_file <- snakemake@output[["counts"]]
dds_file <- snakemake@output[["dds"]]

log_file <- snakemake@log[["log"]]


########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read data
names(count_file_list) <- GenerateName(count_file_list, file = ".htseq-count")
count_data <- rbindlist(lapply(count_file_list, fread), idcol = "library")
setnames(count_data, c("V1", "V2"), c("gene", "reads"))

# generate deseq2 object
count_matrix <- as.matrix(
  data.frame(
    dcast(count_data[!grepl("^__", gene)],
          gene ~ library,
          value.var = "reads"),
    row.names = "gene"))
coldata <- data.frame(row.names = colnames(count_matrix),
           library = colnames(count_matrix))
dds_bg <- DESeqDataSetFromMatrix(count_matrix,
                                 coldata,
                                 design = ~ library)
dds_bg <- estimateSizeFactors(dds_bg)

# write output
fwrite(count_data[!grepl("^__", gene)],
       counts_file)
saveRDS(dds_bg, dds_file)

# write log
sessionInfo()
