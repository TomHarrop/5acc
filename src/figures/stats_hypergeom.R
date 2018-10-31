library(data.table)


tfdb_file <- "output/010_data/tfdb.Rds"
de_genes_file <- "output/050_deseq/wald_tests/expr_genes/all/stage.csv"

alpha <- 0.1

# read files
deseq_res <- fread(de_genes_file)
tfdb <- readRDS(tfdb_file)

# calculate TF enrichment
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

# what about per-family enrichment
families <- tfdb[, c(unique(Family))]
names(families) <- families

all_list <- lapply(families, function(x)
  tfdb[Family == x, unique(`Protein ID`)])
bg_list <- lapply(all_list, function(x)
  data.frame(gene_id = x[x %in% deseq_res$gene_id]))

bg <- rbindlist(bg_list, idcol = "family")

CountFamily <- function(x) {
  deseq_res[gene_id %in% x, length(unique(gene_id))]
}
CountSig <- function(x, alpha) {
  deseq_res[padj < alpha & gene_id %in% x, length(unique(gene_id))]
}

tf_hyper <- bg[, .(n_in_bg = CountFamily(gene_id),
       n_sig = CountSig(gene_id, alpha)),
   by = family]
tf_hyper[, all_in_bg := length(all_in_bg)]
tf_hyper[, all_de := length(all_de)]
tf_hyper[n_sig > 2, 
         p_hyper := phyper(
           q = n_sig - 1,
           m = n_in_bg,
           n = all_in_bg - n_in_bg,
           k = all_de,
           lower.tail = FALSE
         ),
         by = family]
tf_hyper[!is.na(p_hyper), p_adj := p.adjust(p_hyper, "fdr")]
