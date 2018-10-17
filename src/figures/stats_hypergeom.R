library(data.table)


tfdb_file <- "output/010_data/tfdb.Rds"
de_genes_file <- "output/050_deseq/wald_tests/expr_genes/all/stage.csv"
alpha <- 0.1

# read files
deseq_res <- fread(de_genes_file)
tfdb <- readRDS(tfdb_file)

# calculate enrichment
all_tfs <- tfdb[, unique(`Protein ID`)]
all_in_bg <- deseq_res[, unique(gene_id)]
all_de <- deseq_res[padj < alpha, unique(gene_id)]
tfs_de <- intersect(all_tfs, all_de)
tfs_in_bg <- intersect(all_tfs, all_in_bg)

phyper(
  q = length(tfs_de) - 1,
  m = length(tfs_in_bg),
  n = length(all_in_bg) - length(tfs_in_bg),
  k = length(all_de),
  lower.tail = FALSE
)
