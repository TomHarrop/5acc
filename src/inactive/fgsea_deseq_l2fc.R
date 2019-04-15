library(data.table)
library(DESeq2)
library(fgsea)



spec_order <- c("or" = "O. rufipogon",
                "osi" = "O. sativa indica",
                "ob" = "O. barthii",
                "og" = "O. glaberrima")

stage_order <- c(PBM = "IM", SM = "DM")


# data
families_file <- "output/010_data/tfdb_families.Rds"
dds_file <- "output/050_deseq/filtered_dds.Rds"
families <- readRDS(families_file)
dds <- readRDS(dds_file)

###################
# ENRICHMENT TEST #
###################

cpus <- 8
BiocParallel::register(BiocParallel::MulticoreParam(cpus))

# named list of family members
all_fams <- families[, unique(Family)]
names(all_fams) <- all_fams
fam_list <- lapply(all_fams, function(x) 
  families[Family == x, unique(`Protein ID`)])

# try some deseq stats
design(dds) <- ~ accession + stage
dds <- DESeq(dds, parallel = TRUE)
res <- data.table(results(dds,
                          tidy = TRUE,
                          lfcThreshold = log(1.5, 2),
                          alpha = 0.1))

dds_padj <- res[, structure(padj, names = row)]
dds_lfc <- res[, structure(log2FoldChange, names = row)]
dds_stat <- res[, structure(stat, names = row)]

# can we do a general TF enrichment?
fam_list[["tf_list"]] <- unique(unlist(fam_list))

set.seed(1)
fres_padj <- data.table(fgsea(fam_list,
      dds_padj,
      minSize = 15,
      maxSize = length(fam_list$tf_list),
      nperm = 1000,
      nproc = cpus))
setorder(fres_padj, padj)
fres_padj[, 1:7]

set.seed(14)
fres_lfc <- fgsea(fam_list,
      dds_lfc,
      minSize = 15,
      maxSize = length(fam_list$tf_list),
      nperm = 1e5,
      nproc = cpus,
      gseaParam = 2)
fres_lfc_dt <- data.table(fres_lfc)
setorder(fres_lfc_dt, padj)
fres_lfc_dt[, 1:7]

fres_lfc_dt[, unlist(leadingEdge), by = pathway][, .N, by = pathway]

fres_stat <- data.table(fgsea(fam_list,
                              dds_stat,
                              minSize = 15,
                              maxSize = 500,
                              nperm = 1000,
                              nproc = cpus))
setorder(fres_stat, padj)
fres_stat[, 1:7]


plotEnrichment(fam_list$MADS,
               dds_lfc[!is.na(dds_lfc)],
               gseaParam = 1)

plotEnrichment(fam_list$`AP2-EREBP`,
               dds_lfc[!is.na(dds_lfc)],
               gseaParam = 1)


plotGseaTable(fam_list[fres_lfc_dt[1:5, unique(pathway)]],
              dds_lfc[!is.na(dds_lfc)],
              fres_lfc, 
              gseaParam = 1)




