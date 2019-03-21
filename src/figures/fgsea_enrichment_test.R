#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(fgsea)

#############
# FUNCTIONS #
#############

# RunFgsea <- function(pc, pcro, family_list){
#   # vector of scores on pc
#   pc_named <- structure(pcro[[pc]], names = pcro$locus_id)
#   # run fgsea
#   set.seed(14)
#   gsea_res <- fgsea(pathways = fam_list,
#                     stats = pc_named,
#                     minSize = 15,
#                     maxSize = 2000,
#                     nperm = 1000)
#   gsea_dt <- data.table(gsea_res)
#   setorder(gsea_dt, padj)
#   setnames(gsea_dt, "pathway", "family")
#   return(gsea_dt)
# }

###########
# GLOBALS #
###########

families_file <- snakemake@input[["families"]]
wald_stage_file <- snakemake@input[["wald_stage"]]

# plots
fig1_file <- snakemake@output[["fig1"]]
sf1_file <- snakemake@output[["sf1"]]

# tables
tab1_file <- snakemake@output[["table1"]]
tab2_file <- snakemake@output[["table2"]]

cpus <- snakemake@threads[[1]]

# dev
# families_file <- "output/010_data/tfdb_families.Rds"
# wald_stage_file <- "output/050_deseq/wald_tests/expr_genes/all/stage.csv"
# cpus <- 8
# tab1_file <- "tmp/no_pca/fgsea_enrichement.csv"
# tab2_file <- "tmp/no_pca/leading_edge.csv"

spec_order <- c("or" = "O. rufipogon",
                "osi" = "O. sativa indica",
                "ob" = "O. barthii",
                "og" = "O. glaberrima")

stage_order <- c(PBM = "IM", SM = "DM")

########
# MAIN #
########

# data
families <- readRDS(families_file)
wald_stage <- fread(wald_stage_file)


###################
# ENRICHMENT TEST #
###################

# named list of family members
all_fams <- families[, unique(Family)]
names(all_fams) <- all_fams
fam_list <- lapply(all_fams, function(x) 
  families[Family == x, unique(`Protein ID`)])

# get deseq stat to test
deseq_stat <- wald_stage[, structure(abs(log2FoldChange), names = gene_id)]

# subset fam_list, no point running genes that don't exist
fam_list_subset <- lapply(fam_list, function(x) 
  x[x %in% names(deseq_stat)])

# run fgsea
set.seed(1)
fres <- (fgsea(fam_list_subset,
               deseq_stat,
               minSize = 15,
               maxSize = Inf,
               nperm = 1e6,
               nproc = cpus))
fres_dt <- data.table(fres)
setnames(fres_dt, "pathway", "family")
setorder(fres_dt, padj)

# table of genes driving enrichment
leading_edge <- fres_dt[, .(gene_id = unlist(leadingEdge)),
                        by = .(family)]

leading_edge_values <- merge(leading_edge,
                             wald_stage,
                             all.x = TRUE,
                             all.y = FALSE)
setorder(leading_edge_values, padj)

# write output
fwrite(fres_dt[, names(fres_dt) != "leadingEdge", with = FALSE],
       tab1_file)
fwrite(leading_edge_values[, !names(leading_edge_values) %in% c(
  "design", "wald_test", "RapID"),
  with = FALSE],
  tab2_file)

# Log
sessionInfo()



