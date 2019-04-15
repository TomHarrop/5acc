library(data.table)

CountFamily <- function(x, my_res) {
  my_res[gene_id %in% x, length(unique(gene_id))]
}
CountSig <- function(x, my_res, alpha) {
  my_res[padj < alpha & gene_id %in% x, length(unique(gene_id))]
}


tfdb_file <- "output/010_data/tfdb.Rds"
de_genes_file <- "output/050_deseq/wald_tests/expr_genes/all/stage.csv"
# de_genes_file <- "tmp/no_nipponbare/wald/all/stage.csv"
dom_genes_file <- "output/050_deseq/wald_tests/tfs/all/domestication.csv"
# dom_genes_file <- "tmp/no_nipponbare/wald/all/domestication.csv"

alpha <- 0.1

# read files
deseq_res <- fread(de_genes_file)
dom_genes <- fread(dom_genes_file)
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

tf_hyper <- bg[, .(n_in_bg = CountFamily(gene_id, deseq_res),
       n_sig = CountSig(gene_id, deseq_res, alpha)),
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
setorder(tf_hyper, p_adj)

# ap2 enrichment in domestication
dom_bg_list <- lapply(all_list, function(x)
  data.frame(gene_id = x[x %in% dom_genes$gene_id]))

dom_bg <- rbindlist(dom_bg_list, idcol = "family")

dom_all_in_bg = dom_genes[,length(unique(gene_id))]
dom_all_de = dom_genes[padj < alpha,length(unique(gene_id))]


dom_hyper <- dom_bg[, .(n_in_bg = CountFamily(gene_id, dom_genes),
                   n_sig = CountSig(gene_id, dom_genes, alpha)),
               by = family]
dom_hyper[, all_in_bg := dom_all_in_bg]
dom_hyper[, all_de := dom_all_de]
dom_hyper[n_sig > 2, 
         p_hyper := phyper(
           q = n_sig - 1,
           m = n_in_bg,
           n = all_in_bg - n_in_bg,
           k = all_de,
           lower.tail = FALSE
         ),
         by = family]
dom_hyper[!is.na(p_hyper), p_adj := p.adjust(p_hyper, "fdr")]
setorder(dom_hyper, p_adj)

# domestication ap2s in de genes
ap2 <- tfdb[Family == "AP2-EREBP", unique(`Protein ID`)]
ap2_dom <- dom_genes[padj < 0.1 & gene_id %in% ap2, unique(gene_id)]
intersect(ap2_dom, plot_ap2)
deseq_res[gene_id %in% ap2_dom, log2FoldChange]
