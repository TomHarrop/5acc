#!/usr/bin/Rscript

library(data.table)

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

# write output
fwrite(count_data[!grepl("^__", gene)],
       counts_file)

# write log
sessionInfo()
