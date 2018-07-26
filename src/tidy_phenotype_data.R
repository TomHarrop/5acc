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
cali_file <- snakemake@input[["cali"]]


tidy_mtp_file <- snakemake@output[["mtp"]]
tidy_cali_file <- snakemake@output[["cali"]]

# dev
# mtp_file <- "data/phenotyping/Phenotype_PanicleSequenced_corrected2.csv"
# cali_file <- "data/phenotyping/OsOgObOrPTRAPdata_PaperTom.txt"

########
# MAIN #
########


# load the data
pheno_cali <- fread(cali_file, dec = ",")
pheno_mtp <- fread(mtp_file)

# mung cali
cali_species_names <-c(Ob = "barthii",
                       Os = "sativa",
                       Og = "glaberrima",
                       Or = "rufipogon")
cali_vars <- c("RL",
               "PbN",
               "PbL",
               "PbIntL",
               "SbN",
               "SbL",
               "SbIntL",
               "TbN",
               "SpN")
pheno_cali[, Species := factor(plyr::revalue(Origin, cali_species_names),
                               levels = cali_species_names)]

# this can be used for prcomp directly (plot the loadings separately)
tidy_cali <- pheno_cali[, c(
  "Species", "Sowing_nb", "Repet_nb", "Plant_nb", "Panicle_nb",
  cali_vars
), with = FALSE]

# mung mtp
mtp_accession_names <- c(B88 = "barthii",
                         IR64 = "indica",
                         Niponbarre = "japonica",
                         Tog5681 = "glaberrima",
                         W1654 = "rufipogon")
mtp_vars <- c(`Pb_nb (PbN)` = "pbn",
              `Sb_nb (SbN)` = "sbn",
              `Sp_nb (SpN)` = "spn")

names(mtp_accession_names) <- tolower(names(mtp_accession_names))
pheno_mtp[, species_simple := plyr::revalue(tolower(`Accession Name`),
                                           mtp_accession_names)]

setnames(pheno_mtp, names(mtp_vars), mtp_vars)

tidy_mtp <- pheno_mtp[, c("file_name", "species_simple", mtp_vars),
                      with = FALSE]

# write output
fwrite(tidy_cali, tidy_cali_file)
fwrite(tidy_mtp, tidy_mtp_file)

# write log
sessionInfo()

