#!/usr/bin/env Rscript

library(data.table)


#############
# FUNCTIONS #
#############

StripPercent <- function(x){
  as.numeric(sub("%", "", x))
}

DateConvert <- function(x){
  as.POSIXct(strptime(x, format = '%b %d %H:%M:%S'))
}

GenerateName <- function(x, file) {
  paste(
    gsub("^.+/", "", dirname(x)),
    sub(file, "", basename(x)),
    sep = "_")
}


###########
# GLOBALS #
###########

log_file_list <- snakemake@input[["log_file_list"]]
read_count_list <- snakemake@input[["read_count_list"]]

log_file <- snakemake@log[["log"]]

rds_file <- snakemake@output[["rds"]]
csv_file <- snakemake@output[["csv"]]

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read data
names(log_file_list) <- GenerateName(log_file_list, file = ".Log.final.out")
names(read_count_list) <- GenerateName(read_count_list, file = ".ReadsPerGene.out.tab")

star_log_frames <- lapply(log_file_list,
                          read.delim,
                          header = FALSE,
                          sep = "|",
                          fill = TRUE,
                          strip.white = TRUE,
                          stringsAsFactors = FALSE)

read_count_dts <- lapply(read_count_list, fread)

# make into a data.table
rem_cols <- c("UNIQUE READS:",
              "MULTI-MAPPING READS:",
              "UNMAPPED READS:",
              "CHIMERIC READS:")

star_log_data <- rbindlist(star_log_frames, idcol = 'library')
setnames(star_log_data, c("V1", "V2"), c("variable", "value"))

# go wide
star_log_wide <- dcast(star_log_data[!variable %in% rem_cols],
                       library ~ variable)
setnames(star_log_wide,
         old = c("Deletion rate per base", "Insertion rate per base"),
         new = c("Deletion rate per base, %", "Insertion rate per base, %"))

# tidy up columns
num_cols <- c("% of chimeric reads",
              "% of reads mapped to multiple loci",
              "% of reads mapped to too many loci",
              "% of reads unmapped: other",
              "% of reads unmapped: too many mismatches",
              "% of reads unmapped: too short",
              "Average input read length",
              "Average mapped length",
              "Deletion average length",
              "Deletion rate per base, %",
              "Insertion average length",
              "Insertion rate per base, %",
              "Mapping speed, Million of reads per hour",
              "Mismatch rate per base, %",
              "Number of chimeric reads",
              "Number of input reads",
              "Number of reads mapped to multiple loci",
              "Number of reads mapped to too many loci",
              "Number of splices: AT/AC",
              "Number of splices: Annotated (sjdb)",
              "Number of splices: GC/AG",
              "Number of splices: GT/AG",
              "Number of splices: Non-canonical",
              "Number of splices: Total",
              "Uniquely mapped reads %",
              "Uniquely mapped reads number")

date_cols <- c("Started job on",
               "Started mapping on",
               "Finished on")

star_logs_reclassed <- cbind(star_log_wide[, .(library)],
                   star_log_wide[, lapply(.SD, DateConvert),
                                 .SDcols = date_cols],
                   star_log_wide[,lapply(.SD, StripPercent),
                                 .SDcols = num_cols])

# count the number of reads in genes
read_counts <- rbindlist(read_count_dts, idcol = "library")
setnames(read_counts, c("V1", "V4"), c("msuId", "reads"))
read_counts_per_library <- read_counts[!grepl("^N_", msuId),
            .("Number of reads in genes" = sum(reads)),
            by = library]


# merge results
star_logs <- merge(star_logs_reclassed,
                   read_counts_per_library,
                   by = "library")
setkey(star_logs, library)

# save results
saveRDS(star_logs, rds_file)
fwrite(star_logs, csv_file)

# write log
sessionInfo()

