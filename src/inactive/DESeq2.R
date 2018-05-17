#!/usr/bin/Rscript

library(data.table)

# load the quant files
quant.files.all <- list.files("output", pattern = "ReadsPerGene",
                              recursive = TRUE, full.names = TRUE)
quantFiles <- quant.files.all[!grepl("remap", quant.files.all)]
if (length(quantFiles) == 0) {
  stop("Couldn't find quant files, exiting\n")
  quit(save = "no", status = 1)
}
rawCountData <- lapply(quantFiles, read.table, stringsAsFactors = FALSE)

# make a table of countData for DESeq2
names(rawCountData) <- gsub(".*([[:upper:]][[:digit:]]).*", "\\1", quantFiles)
tidyCountData <- function(x) {
  cleanCounts <- data.table::copy(rawCountData[[x]])
  rownames(cleanCounts) <- cleanCounts$V1
  data.table::setnames(cleanCounts, old = "V4", new = x)
  return(subset(cleanCounts, select = x))
}
countDataList <- lapply(names(rawCountData), tidyCountData)
countData <- do.call(cbind, countDataList)

# set up colData
colData.table <- data.table(rn = colnames(countData))
colData.table[, accession := plyr::mapvalues(substr(rn, 1, 1),
                              from = c("J", "I", "R", "G", "B"),
                              to = c("japonica", "indica", "rufipogon",
                                     "glaberrima", "barthii"))]
colData.table[as.numeric(substr(rn, 2, 2)) < 4, stage := "PBM"]
colData.table[as.numeric(substr(rn, 2, 2)) > 3, stage := "SM"]
colData.table[accession %in% c("barthii", "rufipogon"),
              domestication := "wild"]
colData.table[!accession %in% c("barthii", "rufipogon"),
              domestication := "domesticated"]
colData.table[accession %in% c("barthii", "glaberrima"),
              continent := "Africa"]
colData.table[accession %in% c("rufipogon", "indica", "japonica"),
              continent := "Asia"]

colData <- data.frame(colData.table, row.names = "rn")
ind <- sapply(colData, is.character)
colData[ind] <- lapply(colData[ind], factor)
colData$domestication <- relevel(colData$domestication, "wild")
colData$accession <- relevel(colData$accession, "japonica")

# build DESeq2 object for basic comparisons
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = subset(countData, !grepl("^N", rownames(countData))),
  colData = colData,
  design = ~ accession + stage)

# run DESeq2
dds <- DESeq2::DESeq(dds)

# run transformations (need to use blind=TRUE for QC)
vst <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
rld <- DESeq2::rlogTransformation(dds, blind = FALSE)

# normalized counts matrix
normCounts <- DESeq2::counts(dds, normalized = TRUE)

# MAKE OUTPUT FOLDER
outDir <- "output/deseq2"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# save output
saveRDS(dds, paste0(outDir, "/dds.Rds"))
saveRDS(vst, paste0(outDir, "/vst.Rds"))
saveRDS(rld, paste0(outDir, "/rld.Rds"))
saveRDS(normCounts, paste0(outDir, "/normCounts.Rds"))

# SAVE LOGS
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
writeLines(sInf, paste0(outDir, "/SessionInfo.txt"))

# exit 0
quit(save = "no", status = 0)
