#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[["log"]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

###########
# GLOBALS #
###########

cali_file <- snakemake@input[["cali"]]
clusters_file <- snakemake@input[["clusters"]]

correlation_file <- snakemake@output[["correlation"]]
pca_file <- snakemake@output[["pca"]]

# dev
# cali_file <- "test/cali.csv"
# clusters_file <- "output/070_clustering/tfs/annotated_clusters_scaled_l2fc.csv"

########
# MAIN #
########

# load the data
pheno_cali <- fread(cali_file)
clusters <- fread(clusters_file)

# get the mean expression by cluster
accession_order <- c("rufipogon", "indica", "barthii",  "glaberrima")

long_clusters <- melt(clusters,
                      id.vars = c("MsuID","cluster"),
                      measure.vars = accession_order,
                      variable.name = "accession")
long_cores <- long_clusters[, .(core_expression = mean(value)),
                            by = .(accession, cluster)]
cores <- dcast(long_cores, accession ~ cluster)
all_clusters <- long_cores[, as.character(unique(cluster))]

# run PCA on the cali results
cali_vars <- c("rachis_length",
               "primary_branch_number",
               "primary_branch_length",
               "primary_branch_internode_length",
               "secondary_branch_number",
               "secondary_branch_length",
               "secondary_branch_internode_length",
               "tertiary_branch_number",
               "spikelet_number")
pca_data <- pheno_cali[, lapply(.SD, as.numeric), .SDcols = cali_vars]
pc <- prcomp(pca_data, center = TRUE, scale = TRUE)

# bind the PC results and get the mean pc by species
pheno_pc <- cbind(pheno_cali, data.table(pc$x))
pheno_pc_means <- pheno_pc[, lapply(.SD, mean),
                           by = Species,
                           .SDcols = colnames(pc$x)]

# merge the cluster results
cluster_pc <- merge(cores,
                    pheno_pc_means,
                    by.x = "accession",
                    by.y = "Species",
                    all.x = TRUE,
                    all.y = FALSE)

# run the correlations
correlations_wide <- cor(data.frame(cluster_pc,
                                    row.names = "accession"), method = "pearson")
correlations_all <- data.table(melt(correlations_wide,
                                    value.name = "pearson_correlation"))
correlations <- correlations_all[
  Var1 %in% paste0("X", all_clusters) &
    Var2 %in% colnames(pc$x),
  .(cluster = plyr::mapvalues(Var1, paste0("X", all_clusters), all_clusters),
    PC = Var2,
    pearson_correlation)]

# save output
saveRDS(pc, pca_file)
fwrite(correlations, correlation_file)

# write log
sessionInfo()
