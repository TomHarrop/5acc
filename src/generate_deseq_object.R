#!/usr/bin/env Rscript

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

read_count_list <- snakemake@input[["read_count_list"]]

log_file <- snakemake@log[["log"]]

dds_file <- snakemake@output[["dds"]]
vst_file <- snakemake@output[["vst"]]
rld_file <- snakemake@output[["rld"]]
norm_counts_file <- snakemake@output[["norm_counts"]]

# debug
# print(read_count_list)

# dev
# read_count_list <- list.files(path = 'output/030_mapping/star-pass2',
#                               pattern = ".ReadsPerGene.out.tab",
#                               full.names = TRUE,
#                               recursive = TRUE)

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read counts
names(read_count_list) <- GenerateName(read_count_list,
                                       file = ".ReadsPerGene.out.tab")
read_count_dts <- lapply(read_count_list, fread)
read_counts <- rbindlist(read_count_dts, idcol = "library")

# generate count matrix
count_table <- dcast(read_counts[!grepl("^N_", V1)],
                     V1 ~ library,
                     value.var = "V4")
count_matrix <- as.matrix(data.frame(count_table, row.names = "V1"))

# generate sample data
full_names <- c(osj = "japonica",
                osi = "indica",
                or = "rufipogon",
                og = "glaberrima",
                ob = "barthii")

cd_table <- data.table(id = names(read_count_list))
cd_table[, c("spec_code", "stage", "replicate") := tstrsplit(id, "_")]
cd_table[, accession := factor(plyr::revalue(spec_code, full_names),
                               levels = full_names)]
cd_table[accession %in% c("barthii", "rufipogon"),
         domestication := "wild"]
cd_table[!accession %in% c("barthii", "rufipogon"),
         domestication := "domesticated"]
cd_table[, domestication := factor(domestication,
                                   levels = c("wild", "domesticated"))]
cd_table[accession %in% c("barthii", "glaberrima"),
         continent := "Africa"]
cd_table[accession %in% c("rufipogon", "indica", "japonica"),
         continent := "Asia"]

cd <- data.frame(cd_table[, .(id, accession, stage, continent, domestication)],
                 row.names = "id")

# generate deseq object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = cd[colnames(count_matrix),],
                              design = ~ accession + stage)

# run DESeq2
dds <- DESeq2::DESeq(dds,
                     parallel = TRUE)

# run transformations (n.b. use blind=TRUE for QC)
vst <- varianceStabilizingTransformation(dds, blind = FALSE)
rld <- rlogTransformation(dds, blind = FALSE)

# normalized counts matrix
norm_counts <- counts(dds, normalized = TRUE)

# save output
saveRDS(dds, dds_file)
saveRDS(vst, vst_file)
saveRDS(rld, rld_file)
saveRDS(norm_counts, norm_counts_file)

# write log
sessionInfo()
