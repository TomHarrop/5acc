#!/usr/bin/env Rscript

# Table S8: clustered-genes

cluster_file <- "output/070_clustering/tfs/annotated_clusters_scaled_l2fc.csv"

# read data
clusters <- fread(cluster_file)

clusters[, c("RapID",
             "names",
             "OgroObjective",
             "OgroRef") := NULL]
setnames(clusters, "MsuID", "gene_id")
setcolorder(clusters, c("cluster", "gene_id"))

fwrite(clusters, "test/Table_S8.csv")
