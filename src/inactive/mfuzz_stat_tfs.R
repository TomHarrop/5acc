#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[["log"]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


library(data.table)
library(DESeq2)
library(ggplot2)
library(Mfuzz)

#############
# FUNCTIONS #
#############

CountGenesOfInterest <- function(ids, genes_of_interest){
  sum(unique(ids) %in% genes_of_interest)
}

ClusterCounter <- function(x) {
  cluster_table[cluster == x, length(unique(MsuID))]
}

ExtractAccessionResults <- function(accession, dds, lfcThreshold, alpha) {
  my_num <- paste(accession, "SM", sep = "_")
  my_den <- paste(accession, "PBM", sep = "_")
  my_res <- results(dds,
                    contrast = c("group", my_num, my_den),
                    lfcThreshold = lfcThreshold,
                    alpha = alpha)
  my_dt <- data.table(data.frame(my_res), keep.rownames = TRUE)
  my_dt[, .(gene_id = rn, stat)]
}


###########
# GLOBALS #
###########

tfdb_file <- snakemake@input[["tfdb"]]
dds_file <- snakemake@input[["dds"]]

alpha <- snakemake@params[["alpha"]]
lfcThreshold <- snakemake@params[["lfc_threshold"]]
seed <- snakemake@params[["seed"]]

cpus <- snakemake@threads[[1]]

log_file <- snakemake@log[["log"]]

cluster_plot_file <- snakemake@output[["cluster_plot"]]
hyper_test_file <- snakemake@output[["hyper"]]
cluster_file <- snakemake@output[["clusters"]]

# dev 
# dds_file <- "output/050_deseq/dds_tfs.Rds"
# tfdb_file <- "output/010_data/tfdb.Rds"
# alpha <- 0.1
# lfcThreshold <- log(1.5, 2)
# cpus <- 8
# seed <- 1


########
# MAIN #
########

BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# read the data
tfdb <- readRDS(tfdb_file)
dds <- readRDS(dds_file)

# list of AP2s
ap2_genes <- tfdb[Family == "AP2-EREBP", unique(`Protein ID`)]

# re-run DESeq2 by group
dds$group <- as.factor(paste(dds$accession, dds$stage, sep = "_"))
design(dds) <-  ~ group
dds <- DESeq(dds, parallel = TRUE)

# extract results for DE between stages in each accession
accessions <- c(
#  "japonica" = "japonica",
  "indica" = "indica",
  "rufipogon" = "rufipogon",
  "barthii" = "barthii",
  "glaberrima" = "glaberrima")
acc_results_list <- lapply(accessions,
                           ExtractAccessionResults,
                           dds,
                           lfcThreshold,
                           alpha)
acc_results <- rbindlist(acc_results_list, idcol = "accession")
stat_matrix <- as.matrix(
  data.frame(
    dcast(acc_results, gene_id ~ accession),
    row.names = "gene_id"))

# prepare mfuzz object
pheno_data <- data.frame(row.names = colnames(stat_matrix),
                         stage = colnames(stat_matrix))
vg <- ExpressionSet(assayData = stat_matrix,
                    phenoData = new('AnnotatedDataFrame', data = pheno_data))

# standardise
vg_s <- standardise(vg)
vg_f <- filter.NA(vg_s,thres = 0)

# run the clustering
message(paste("Clustering with seed", seed))
set.seed(seed)
c1 <- mfuzz(vg_f, c = 9, m = 1.6)
clusters <- acore(vg_f, c1, min.acore = 0.7)

# reformat clustering results
cluster_list <- lapply(clusters, function(x) data.table(MsuID = rownames(x)))
cluster_table <- rbindlist(cluster_list, idcol = "cluster")

# plot data
pd_wide <- merge(cluster_table,
                 data.table(as.data.frame(exprs(vg_s)), keep.rownames = TRUE),
                 all.x = TRUE,
                 all.y = FALSE,
                 by.x = "MsuID",
                 by.y = "rn")
pd <- melt(pd_wide,
           id.vars = c("MsuID", "cluster"),
           variable.name = "accession",
           value.name = "Scaled z-statistic")

# labels and order
pd[, n_genes := length(unique(MsuID)), by = cluster]
pd[, c_lab := paste0("Cluster ", cluster, " (", n_genes, " genes)")]
accession_order <- c(
  "rufipogon",
  "indica",
  "barthii",
  "glaberrima")
pd[, accession := factor(accession, levels = accession_order)]
setkey(pd, accession, MsuID)

# plot
gp <- ggplot(pd, aes(x = accession, y = `Scaled z-statistic`, group = MsuID)) +
  theme_minimal() +
  xlab(NULL) +
  facet_wrap(~c_lab, ncol = 3) +
  geom_path()

# annotate
annotation_table <- dcast(pd,
                          cluster + MsuID ~ accession,
                          value.var = "Scaled z-statistic")
annotations <- annotation_table[, oryzr::LocToGeneName(unique(MsuID))]
annotated_clusters <- merge(annotation_table,
                            annotations,
                            by = "MsuID",
                            all.x = TRUE,
                            all.y = FALSE)
setkey(annotated_clusters, cluster, MsuID)

# Hypergeometric test of per-family enrichment
families <- tfdb[, c(unique(Family))]
names(families) <- families

# count the number of members in each cluster
family_counts_wide <- cluster_table[, lapply(families, function(x)
  CountGenesOfInterest(MsuID, tfdb[Family == x, unique(`Protein ID`)])),
  by = cluster]
family_counts <- melt(family_counts_wide,
                      id.vars = "cluster",
                      variable.name = "family",
                      value.name = "n_in_cluster")

# get  cluster size
family_counts[, all_in_cluster := ClusterCounter(cluster), by = cluster]

# background size (all genes submitted for clustering)
family_counts[, all_in_bg := length(unique(rownames(vg_s)))]

# merge in number of members in background
bg_counts <- tfdb[
  , .(n_in_bg = CountGenesOfInterest(
    unique(rownames(vg_s)),
    unique(`Protein ID`))),
  by = Family]
family_hyper <- merge(family_counts,
                      bg_counts,
                      by.x = "family",
                      by.y = "Family")

# run a hypergeometric test only if more than one member was clustered
family_hyper[n_in_cluster > 1, p_hyper := phyper(n_in_cluster - 1,
                                                 n_in_bg,
                                                 all_in_bg - n_in_bg,
                                                 all_in_cluster,
                                                 lower.tail = FALSE)]

# adjust the p-values for multiple trials
family_hyper[!is.na(p_hyper), p_adj := p.adjust(p_hyper, "fdr")]

# sort the table
setorder(family_hyper, cluster, p_adj, family, na.last = TRUE)

# save output
ggsave(cluster_plot_file,
       gp,
       device = cairo_pdf,
       width = 10,
       height = 7.5,
       units = "in")
fwrite(family_hyper, hyper_test_file)
fwrite(annotated_clusters, cluster_file)

# write log
sessionInfo()

