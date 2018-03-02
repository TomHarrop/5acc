#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)
library(ggplot2)

#############
# FUNCTIONS #
#############

AnnotateDESeqResultsTable <- function(results_set) {
  my_results <- data.table(as.data.frame(results_set), keep.rownames = TRUE)
  setnames(my_results, "rn", "MsuID")
  my_annotation <- my_results[, oryzr::LocToGeneName(unique(MsuID))]
  merge(my_results, my_annotation, by = "MsuID", all.x = TRUE)
}


PrintAnnotatedDESeqResultsTable <- function(annotated_results_table,
                                            file_name,
                                            outdir) {
  my_outfile <- file.path(outdir,
                          paste0(file_name, ".tab"))
  fwrite(annotated_results_table,
         my_outfile,
         sep = "\t")  
}  

###########
# GLOBALS #
###########

tfdb_file <- "data/tfdb/tfdb_os.Rds"
det_genes_file <- "output/tpm/detected_genes.Rds"
dds_file <- "output/deseq2/dds.Rds"
cpus <- 8
alpha <- 0.1
lfcThreshold <- log(1.5, 2)
outdir <- "output/deseq2/wald_tests/tf_only"

########
# MAIN #
########

BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# read the dds
dds <- readRDS(dds_file)

# get the list of detected genes
det_genes <- readRDS(det_genes_file)

# get a list of transcription factors
tfdb <- readRDS(tfdb_file)
all_tfs <- tfdb[, unique(Protein.ID)]

# subset the dds for expressed TFs
kept_genes <- intersect(names(dds), intersect(all_tfs, det_genes))
dds_exp_tfs <- dds[kept_genes]

# remove japonica
dds_no_jp <- dds_exp_tfs[, dds_exp_tfs$accession != "japonica"]
dds_no_jp$accession <- droplevels(dds_no_jp$accession)

# DE analysis
dds_basal <- copy(dds_exp_tfs)
dds_stage_continent <- copy(dds_no_jp)
dds_wild_vs_dom_by_group <- copy(dds_no_jp)
dds_domestication <- copy(dds_no_jp)

design(dds_basal) <- ~ accession + stage
dds_basal <- DESeq(dds_basal, parallel = TRUE)

design(dds_stage_continent) <- ~ stage + continent + stage:continent
dds_stage_continent <- DESeq(dds_stage_continent, parallel = TRUE)

design(dds_wild_vs_dom_by_group) <- ~ stage + accession + stage:accession
dds_wild_vs_dom_by_group <- DESeq(dds_wild_vs_dom_by_group, parallel = TRUE)
dds_wild_vs_dom_by_group$accession <- 
  relevel(dds_wild_vs_dom_by_group$accession, "barthii") # barthii must be the base level

design(dds_domestication) <- ~ stage + domestication + stage:domestication
dds_domestication <- DESeq(dds_domestication, parallel = TRUE)

# basal branching TFs
results_basal <- results(dds_basal,
                         contrast = c("stage", "SM", "PBM"),
                         lfcThreshold = lfcThreshold,
                         alpha = alpha)

# de between stages is different for asian and african rice
results_stage_continent <- results(dds_stage_continent,
                                   name = "stageSM.continentAsia",
                                   alpha = alpha)

# accession-specific domestication differences
results_glab_barthii <- results(
  dds_wild_vs_dom_by_group,
  name = "stageSM.accessionglaberrima",
  alpha = alpha)

results_indica_rufipogon <- results(
  dds_wild_vs_dom_by_group,
  contrast = list("stageSM.accessionindica",
                  "stageSM.accessionrufipogon"),
  alpha = alpha)

# domestication / stage interaction
results_domestication <- results(
  dds_domestication,
  name = "stageSM.domesticationdomesticated",
  alpha = alpha)

# collect the results
results_list <- c(
  "fig3_stage" = subset(results_basal, padj < alpha),
  "fig4_interaction_stage_continent" = subset(results_stage_continent, padj < alpha),
  "fig5a_interaction_stage_accession_indica" = subset(results_indica_rufipogon, padj < alpha),
  "fig5b_interaction_stage_accession_glaberrima" = subset(results_glab_barthii, padj < alpha),
  "fig6_domestication" = subset(results_domestication, padj < alpha))

# annotate
annotated_results_list <- lapply(results_list, AnnotateDESeqResultsTable)

# write results
sapply(names(annotated_results_list), function(x)
       PrintAnnotatedDESeqResultsTable(annotated_results_list[[x]],
       file_name = x,
       outdir = outdir))

saveRDS(dds_basal, file.path(outdir, "dds_tf_only.Rds"))

# quick plot
results_domestication[order(results_domestication$padj),]
my_gene <- "LOC_Os05g41780"
oryzr::LocToGeneName(my_gene)
plotCounts(dds_domestication, my_gene, intgroup = c("accession", "stage"))
my_counts <- counts(dds_domestication, normalized = TRUE)[my_gene,]
pd <- data.table(id = my_gene,
           sample = names(my_counts),
           my_counts)
pd[, accession := gsub("[[:digit:]]+", "", sample)]
pd[, sample_number := as.numeric(
  gsub("[^[:digit:]]+", "", sample))]
pd[ sample_number <= 3, stage := "PBM"]
pd[ sample_number >= 4, stage := "SM"]
accession_order <- c("R" = "rufipogon",
                     "I" = "indica",
                     "B" = "barthii",
                     "G" = "glaberrima")
pd[, accession := factor(plyr::revalue(accession, accession_order),
                         levels = accession_order)]
gp <- ggplot(pd, aes(x = stage, y = my_counts, group = accession)) + 
  theme_minimal() +
  xlab("Normalized count") + ylab(NULL) +
  ggtitle(oryzr::LocToGeneName(my_gene)[, unique(symbols)]) +
  facet_wrap(~ accession) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(position = position_jitter(width = 0.2))

ggsave(file.path(outdir, paste0(my_gene, ".pdf")), gp, width = 5, height = 3.75)
