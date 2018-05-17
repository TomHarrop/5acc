library(DESeq2)
library(data.table)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

GenerateMessage(paste("Likelihood ratio test for interaction between stage",
                      "and tertiary axis number"))

# set up parallel processing
cpus <- as.numeric(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
if (is.na(cpus)) {
  cpus <- 1
}
GenerateMessage(paste("Allocating", cpus, "cpu(s)"))
BiocParallel::register(BiocParallel::MulticoreParam(cpus))

design.full <- ~ ta.number + stage + ta.number:stage
design.reduced <- ~ ta.number + stage

# load data
GenerateMessage("Loading base DESeq2 object")
dds.file <- "output/deseq2/dds.Rds"
if (!file.exists(dds.file)) {
  stop("Couldn't find dds.Rds")
}
dds <- readRDS(dds.file)

# pre-filter based on tpm cutoffs
GenerateMessage("Removing undetected genes")
detected.genes.file <- "output/tpm/detected_genes.Rds"
if (!file.exists(detected.genes.file)) {
  stop("Couldn't find detected_genes.Rds")
}
detected.genes <- readRDS(detected.genes.file)
counts.table <- data.table(DESeq2::counts(dds)[detected.genes, ],
                           keep.rownames = TRUE, key = "rn")

# set up ta_number variable
GenerateMessage("Loading phenotype data")
pt.file <- "data/phenotyping/Phenotype_Panicle_corrected.csv"
if (!file.exists(pt.file)) {
  stop("Couldn't find Phenotype_Panicle_corrected.csv")
}
pt.data <- data.table(read.csv(pt.file, stringsAsFactors = FALSE))
GenerateMessage("Setting up numeric ta.number variable")
pt.data <- pt.data[!grep("Echelle", file_name)]
TidyFileName <- function(x) {
  toupper(gsub("^([[:alnum:]]+).*", "\\1", x))
}
pt.data[, accession := plyr::mapvalues(
  TidyFileName(file_name),
  from = c("B88", "IR64", "NIP", "TOG5681", "W1654"),
  to = c("barthii", "indica", "japonica", "glaberrima", "rufipogon"))]
mean.ta.number <- pt.data[, .(ta.number = mean(TA_nb)), by = accession]
#mean.ta.number <- pt.data[, .(ta.number = mean(Sp_nb)), by = accession]
mean.ta.number[, ta.number := cut(ta.number, breaks = 3)]
setkey(mean.ta.number, "accession")

# set up new deseq object
GenerateMessage("Regenerating DESeq2 object")
colData <- data.table(data.frame(SummarizedExperiment::colData(dds)),
                      keep.rownames = TRUE, key = "accession")
col.data.table <- mean.ta.number[colData, .(
  rn, accession, ta.number, stage, sizeFactor)]
setkey(col.data.table, "rn")
col.data <- data.frame(col.data.table, row.names = "rn")
setcolorder(counts.table, c("rn", rownames(col.data)))
counts <- data.frame(counts.table, row.names = "rn")
dds.ta.number <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts, design = design.full, colData = col.data)

# run LRT
GenerateMessage("Running the LRT")
dds.ta.number <- DESeq2::DESeq(dds.ta.number, test = "LRT",
                               reduced = design.reduced, parallel = TRUE)

# explore results
resultsNames(dds.ta.number)
tpm <- readRDS("output/tpm/tpm_with_calls.Rds")
res <- results(dds.ta.number, alpha = 0.05, parallel = TRUE)
summary(res)



sig.genes <- tpm[rownames(res[order(res$padj), ][c(1:20),])]

res.sig <- subset(res, padj < 0.05)
sig.genes <- tpm[rownames(res.sig[order(abs(res.sig$lfcSE), decreasing = TRUE), ][c(1:20),])]

sig.genes[, gene.name := oryzr::LocToGeneName(gene)$symbols, by = gene]
sig.genes[is.na(gene.name), gene.name := gene]
ggplot(sig.genes, aes(x = stage, y = tpm, colour = species, group = species)) +
  facet_wrap(~gene.name, scales = "free_y") +
  scale_colour_brewer(palette = "Set1") +
  geom_smooth(method = "lm", se = FALSE, size = 0.5, alpha = 0.5) + 
  geom_point(position = position_jitter(width = 0.2))



