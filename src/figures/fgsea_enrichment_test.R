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

RunFgsea <- function(pc, pcro, family_list){
  # vector of scores on pc
  pc_named <- structure(pcro[[pc]], names = pcro$locus_id)
  # run fgsea
  gsea_res <- fgsea(pathways = fam_list,
                    stats = pc_named,
                    minSize = 15,
                    maxSize = 500,
                    nperm = 1000)
  gsea_dt <- data.table(gsea_res)
  setorder(gsea_dt, padj)
  setnames(gsea_dt, "pathway", "family")
  return(gsea_dt)
}

###########
# GLOBALS #
###########

families_file <- snakemake@input[["families"]]
pcro_file <- snakemake@input[["pcro"]]

# plots
fig1_file <- snakemake@output[["fig1"]]
sf1_file <- snakemake@output[["sf1"]]

# tables
tab1_file <- snakemake@output[["table1"]]
tab2_file <- snakemake@output[["table2"]]
tab3_file <- snakemake@output[["table3"]]


# dev
# families_file <- "output/010_data/tfdb_families.Rds"
# pcro_file <- "output/050_deseq/rlog_pca/pcro.Rds"
# pcro_file <- "tmp/no_nipponbare/rlog_pca/pcro.Rds"
# tab1_file <- "tmp/no_nipponbare/rlog_pca/pca_enrichment.csv"
# tab2_file <- "tmp/no_nipponbare/rlog_pca/genes_driving_enrichment.csv"
# tab3_file <- "tmp/no_nipponbare/rlog_pca/pca_rotation.csv"

spec_order <- c("or" = "O. rufipogon",
                "osi" = "O. sativa indica",
                "ob" = "O. barthii",
                "og" = "O. glaberrima")

stage_order <- c(PBM = "IM", SM = "DM")

########
# MAIN #
########

# data
pcro <- data.table(readRDS(pcro_file))
families <- readRDS(families_file)

###################
# ENRICHMENT TEST #
###################

pcs_of_interest <- paste0("PC", 1:4)
names(pcs_of_interest) <- pcs_of_interest

# named list of family members
all_fams <- families[, unique(Family)]
names(all_fams) <- all_fams
fam_list <- lapply(all_fams, function(x) 
  families[Family == x, unique(`Protein ID`)])

# tidy pcro table
pcro_long <- melt(pcro, id.vars = "locus_id",
                  measure.vars = pcs_of_interest,
                  variable.name = "component")
pcro_long[, abs_val := abs(value)]
setorder(pcro_long, component, -abs_val)
pcro_long[, abs_rank := 1:length(unique(locus_id)), by = component]
setorder(pcro_long, component, abs_rank)

# run fgsea on all PCS of interest
gsea_list <- lapply(pcs_of_interest,
                    RunFgsea,
                    pcro = pcro,
                    family_list = fam_list)
gsea_all <- rbindlist(gsea_list, idcol = "component")
setorder(gsea_all, component, padj)

# table of genes driving enrichment
leading_edge <- gsea_all[, .(locus_id = unlist(leadingEdge)),
                         by = .(family, component)]

leading_edge_values <- merge(leading_edge,
                             pcro_long,
                             all.x = TRUE,
                             all.y = FALSE)

# add annotations
annot <- pcro[, oryzr::LocToGeneName(unique(locus_id))]
leading_edge_annot <- merge(leading_edge_values,
                            annot[, .(locus_id = MsuID, symbols, names, MsuAnnotation)],
                            by = "locus_id",
                            all.x = TRUE,
                            all.y = FALSE)
setorder(leading_edge_annot, component, abs_rank)

pcro_annot <- merge(pcro,
      annot[, .(locus_id = MsuID, symbols, names, MsuAnnotation)],
      by = "locus_id",
      all.x = TRUE,
      all.y = FALSE)

# write output
fwrite(gsea_all[, names(gsea_all) != "leadingEdge", with = FALSE],
       tab1_file)
fwrite(leading_edge_annot[, names(leading_edge_annot) != "abs_val",
                          with = FALSE],
       tab2_file)
fwrite(pcro_annot, tab3_file)

# Log
sessionInfo()



