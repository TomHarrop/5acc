#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)
library(Mfuzz)

#############
# FUNCTIONS #
#############

CountGenesOfInterest <- function(ids, genes_of_interest){
  sum(unique(ids) %in% genes_of_interest)
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
set.seed(1)


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
  #"japonica" = "japonica",
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
c1 <- mfuzz(vg_s, c = 6, m = 2.1)
clusters <- acore(vg_s, c1, min.acore = 0.7)
mfuzz.plot(vg_s, cl = c1, mfrow = c(2, 3), time.labels = colnames(vg_s), min.mem = 0.7)

sapply(clusters, function(x) length(unique(rownames(x))))
sapply(clusters, function(x)
  CountGenesOfInterest(unique(rownames(x)), ap2_genes))
oryzr::LocToGeneName(unique(rownames(clusters[[3]])))
oryzr::LocToGeneName(unique(rownames(clusters[[3]])))[, unique(symbols)]
