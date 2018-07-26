#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[["log"]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

###########
# GLOBALS #
###########

mtp_file <- snakemake@input[["mtp"]]
clusters_file <- snakemake@input[["clusters"]]

correlation_file <- snakemake@output[["correlation"]]

# dev
# mtp_file <- "output/080_phenotype/mtp.csv"
# clusters_file <- "output/070_clustering/tfs/annotated_clusters_scaled_l2fc.csv"

########
# MAIN #
########

# load the data
pheno_mtp <- fread(mtp_file)
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

# can we do a PCA on this... no, not enough indivs
pc <- pheno_mtp[, prcomp(
  data.table(pbn, sbn, spn),
  center = TRUE, scale = TRUE
)]
pc_results <- cbind(pheno_mtp, pc$x)

# summarise the phenotype results
pheno_vars <- c("pbn", "sbn", "spn")
median_pt <- pheno_mtp[, lapply(.SD, median, na.rm = TRUE),
                       by = species_simple, .SDcols = pheno_vars]

# join
pheno_cores <- merge(cores,
                     median_pt,
                     by.x = "accession",
                     by.y = "species_simple")

# calculate correlations
pbn_corr <- pheno_cores[, lapply(.SD, function(x)
  cor(x, pbn, method = "pearson")),
  .SDcols = all_clusters]
sbn_corr <- pheno_cores[, lapply(.SD, function(x)
  cor(x, sbn, method = "pearson")),
  .SDcols = all_clusters]
spn_corr <- pheno_cores[, lapply(.SD, function(x)
  cor(x, spn, method = "pearson")),
  .SDcols = all_clusters]

all_corr_wide <- rbindlist(list(pbn = pbn_corr,
                                sbn = sbn_corr,
                                spn = spn_corr),
                           idcol = "variable")
all_corr <- melt(all_corr_wide,
                 id.vars = "variable",
                 variable.name = "cluster",
                 value.name = "pearson_correlation")

# save output
fwrite(all_corr, correlation_file)

# write log
sessionInfo()
