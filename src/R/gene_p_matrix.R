#!/usr/bin/env Rscript

library(data.table)

# load wald test results
wald.files <- list.files("output/deseq2/wald_tests",
                         pattern = "results_table.Rds", full.names = TRUE,
                         include.dirs = FALSE)
names(wald.files) <- gsub("_results_table.Rds", "",
                          basename(wald.files), fixed = TRUE)
wald.tables <- lapply(wald.files, readRDS)

# make gene to p list (long)
wald.results <- rbindlist(wald.tables, idcol = "wald_test", fill = TRUE)
wald.results[!is.na(domestication),
             wald_test := paste(wald_test, domestication, sep = "_")]
wald.results[!is.na(continent),
             wald_test := paste(wald_test, continent, sep = "_")]
wald.results[!is.na(accession),
             wald_test := paste(wald_test, accession, sep = "_")]
wald.results[!is.na(contrast),
             wald_test := paste(wald_test, contrast, sep = "_")]

# make gene vs. p table
gene.p.table <- dcast.data.table(wald.results,
                                 gene ~ wald_test,
                                 value.var = "padj")

# find genes in QTLs and add them to table
gwas.regions <- readRDS(
  "/home/tom/Desktop/restored/Projects/gwas-combine/output/regions_with_genes.Rds")
crowell.genes <- gwas.regions[!is.na(crowell.trait) & crowell.trait != "", unique(gene)]
lmi.genes <- gwas.regions[!is.na(lmi.trait) & lmi.trait != "", unique(gene)]

gene.p.table[gene %in% crowell.genes, crowell_qtl := TRUE]
gene.p.table[is.na(crowell_qtl), crowell_qtl := FALSE]

gene.p.table[gene %in% lmi.genes, lmi_qtl := TRUE]
gene.p.table[is.na(lmi_qtl), lmi_qtl := FALSE]

# add annotation (slow way)
annotation <- gene.p.table[, oryzr::LocToGeneName(unique(gene))]
annotated.p.table <- merge(gene.p.table,
                           annotation[, .(MsuID, symbols)],
                           by.x = "gene",
                           by.y = "MsuID",
                           all.x = TRUE)
c.order <- c("gene",
             "symbols",
             "domestication",
             "domestication_asia_indica",
             "domestication_asia_japonica",
             "domestication_by_continent_africa",
             "domestication_by_continent_asia",
             "species_barthii.glaberrima",
             "species_barthii.indica",
             "species_barthii.rufipogon",
             "species_glaberrima.indica",
             "species_glaberrima.rufipogon",
             "species_indica.rufipogon",
             "species_japonica.barthii",
             "species_japonica.glaberrima",
             "species_japonica.indica",
             "species_japonica.rufipogon",
             "stage",
             "stage_continent_africa",
             "stage_continent_asia",
             "stage_species_barthii",
             "stage_species_glaberrima",
             "stage_species_indica",
             "stage_species_japonica",
             "stage_species_rufipogon",
             "crowell_qtl",
             "lmi_qtl")
setcolorder(annotated.p.table, c.order)

###############
# Save output #
###############

rutils::GenerateMessage("Writing output")

outdir = "output/gene_p_matrix"
if(!dir.exists(outdir)){
  dir.create(outdir)
}
saveRDS(annotated.p.table, paste0(outdir, "/gene_p_matrix.RDS"))
write.table(annotated.p.table,
            paste0(outdir, "/gene_p_matrix.csv"),
            sep = ",",
            quote = FALSE,
            na = "",
            row.names = FALSE)

sinfo <- rutils::GitSessionInfo()
writeLines(sinfo, paste0(outdir, "/SessionInfo.gene_p_matrix.txt"))