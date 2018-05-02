#!/usr/bin/env Rscript

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
  my_dt[, .(gene_id = rn, log2FoldChange)]
}

MinUnclusteredAp2s <- function(es, m, c, memCutoff) {
  # run the clustering
  my_flclust <- mfuzz(es, c = c, m = m)
  my_acore <- acore(es, my_flclust, min.acore = memCutoff)
  # count the number of clustered genes
  my_sums <- sapply(my_acore, function(x) length(unique(rownames(x))))
  cluster_score <- sum(my_sums) / length(unique(rownames(es)))
  print(my_sums)
  print(cluster_score)
  # count the number of ap2 genes
  my_numbers <- sapply(my_acore, function(x) 
    CountGenesOfInterest(unique(rownames(x)),
                         genes_of_interest = ap2_genes))
  ap2_score <- max(my_numbers / sum(my_numbers))
  print(my_numbers)
  print(ap2_score)
  # return weighted score
  (2 * ap2_score) + cluster_score
  
}

CountClusteredAp2s <- function(es, m, c, memCutoff){
  # run the clustering
  my_flclust <- mfuzz(es, c = c, m = m)
  my_acore <- acore(es, my_flclust, min.acore = memCutoff)
  # return the stats
  data.table(
    clustered_genes = sapply(my_acore, function(x) length(unique(rownames(x)))),
    ap2_genes = sapply(my_acore, function(x) 
      CountGenesOfInterest(unique(rownames(x)),
                           genes_of_interest = ap2_genes)),
    memCutoff = memCutoff)
}


###########
# GLOBALS #
###########

dds_file <- "output/deseq2/wald_tests/tf_only/dds_tf_only.Rds"
tfdb_file <- "data/tfdb/tfdb_os.Rds"
alpha <- 0.1
lfcThreshold <- log(1.5, 2)
cpus <- 8
seed <- 42
outdir <- "output/deseq2/wald_tests/tf_only" # ! GROSS !


########
# MAIN #
########

BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# read the tfdb
tfdb <- readRDS(tfdb_file)
ap2_genes <- tfdb[Family == "AP2-EREBP", unique(Protein.ID)]

dds <- readRDS(dds_file)

# re-run DESeq2 by group
dds$group <- as.factor(paste(dds$accession, dds$stage, sep = "_"))
design(dds) <-  ~ group
dds <- DESeq(dds, parallel = TRUE)

# extract results
accessions <- c(
  # "japonica" = "japonica",
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
lfc_matrix <- as.matrix(
  data.frame(
    dcast(acc_results, gene_id ~ accession),
    row.names = "gene_id"))

# prepare mfuzz object
pd <- data.frame(row.names = colnames(lfc_matrix),
                 stage = colnames(lfc_matrix))
vg <- ExpressionSet(assayData = lfc_matrix,
                    phenoData = new('AnnotatedDataFrame', data = pd))

# standardise
vg_s <- standardise(vg)

# optimise m/c
trial_table <- data.table(
  expand.grid(c(4:15),
              seq(1.5, 2.5, 0.1)))[, .(c = Var1, m = Var2)]
trial_results <- trial_table[, CountClusteredAp2s(vg_s, m, c, 0.5),
                             by = .(m, c)]

trial_results[, prop_ap2 := ap2_genes / clustered_genes]

pd <- trial_results[, .(
  total_clustered_genes = sum(clustered_genes),
  max_prop_ap2 = max(prop_ap2)), by = .(m, c)]
ggplot(pd, aes(x = total_clustered_genes, y = max_prop_ap2)) +
  geom_smooth() +
  geom_point()
pd[max_prop_ap2 > 0.16][which.max(total_clustered_genes)]

# run the clustering
set.seed(seed)
c1 <- mfuzz(vg_s, c = 6, m = 2.1)
clusters <- acore(vg_s, c1, min.acore = 0.7)

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
           value.name = "Scaled L2FC")

# labels and order
pd[, n_genes := length(unique(MsuID)), by = cluster]
pd[, c_lab := paste0("Cluster ", cluster, " (", n_genes, " genes)")]

accession_order <- c(
  "rufipogon",
  "indica",
  # "japonica",
  "barthii",
  "glaberrima")
pd[, accession := factor(accession, levels = accession_order)]

setkey(pd, accession, MsuID)

# plot
gp <- ggplot(pd, aes(x = accession, y = `Scaled L2FC`, group = MsuID)) +
  theme_minimal() +
  xlab(NULL) +
  facet_wrap(~c_lab, ncol = 3) +
  geom_path()

# annotate
annotations <- cluster_table[, oryzr::LocToGeneName(unique(MsuID))]
annotated_clusters <- merge(cluster_table,
                            annotations,
                            by = "MsuID",
                            all.x = TRUE,
                            all.y = FALSE)
setkey(annotated_clusters, cluster, MsuID)

# counts of ap2
enrichment <- cluster_table[, .(
  n_genes = length(unique(MsuID)),
  n_ap2 = CountGenesOfInterest(unique(MsuID), ap2_genes),
  all_genes = length(unique(rownames(vg_s))),
  all_ap2s = CountGenesOfInterest(unique(rownames(vg_s)), ap2_genes)
), by = cluster]
enrichment[, p_hyper := phyper(n_ap2 - 1,
                               all_ap2s,
                               all_genes - all_ap2s,
                               n_genes,
                               lower.tail = FALSE)]
enrichment[, p_adj := p.adjust(p_hyper, "fdr")]

# counts of all families
families <- tfdb[, c(unique(Family))]
names(families) <- families

family_counts_wide <- cluster_table[, lapply(families, function(x)
  CountGenesOfInterest(MsuID, tfdb[Family == x, unique(Protein.ID)])),
  by = cluster]
family_counts <- melt(family_counts_wide,
                      id.vars = "cluster",
                      variable.name = "family",
                      value.name = "n_in_cluster")

family_counts[, all_in_cluster := ClusterCounter(cluster), by = cluster]
family_counts[, all_in_bg := length(unique(rownames(vg_s)))]

bg_counts <- tfdb[, .(n_in_bg = CountGenesOfInterest(unique(rownames(vg_s)),
                                                  unique(Protein.ID))),
                  by = Family]
family_hyper <- merge(family_counts, bg_counts, by.x = "family", by.y = "Family")

family_hyper[n_in_cluster > 1, p_hyper := phyper(n_in_cluster - 1,
                                 n_in_bg,
                                 all_in_bg - n_in_bg,
                                 all_in_cluster,
                               lower.tail = FALSE)]
family_hyper[!is.na(p_hyper), p_adj := p.adjust(p_hyper, "fdr")]
setorder(family_hyper, cluster, family)
family_hyper[p_adj < 0.05]

# write output
fwrite(annotated_clusters,
       file.path(outdir, "clusters.tab"),
       sep = "\t")
fwrite(enrichment,
       file.path(outdir, "cluster_ap2_hypergeometric_test.csv"))
fwrite(family_hyper,
       file.path(outdir, "cluster_all-families_hypergeometric_test.csv"))
ggsave(file.path(outdir, "clusters.pdf"),
       gp, width = 10, height = 7.5)

