# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

###########
# GLOBALS #
###########

raw_tfdb_file <- snakemake@input[["tfdb"]]
raw_families_file <- snakemake@input[["families"]]
output_file <- snakemake@input[["tfdb_with_families"]]

# dev
# raw_tfdb_file <- "data/genome/os/tfdb.tab"
# raw_families_file <- 'data/genome/os/Osj_TF_list'

########
# MAIN #
########

# read data
tfdb_raw <- fread(raw_tfdb_file)
families_raw <- fread(raw_families_file)

# remove duplicates
unique_tfs <- unique(tfdb_raw, by = c("Species", "Family", "Protein ID"))

# extract OS TFs
os <- unique(unique_tfs[Species == 'Oryza sativa subsp. japonica', .(
    "Protein ID" = gsub("\\..*", '', `Protein ID`),
    Family
)])

# manually add ALOG genes
ALOG <- c('LOC_Os07g04670', 'LOC_Os02g07030', 'LOC_Os06g46030',
          'LOC_Os02g41460', 'LOC_Os04g43580', 'LOC_Os10g33780',
          'LOC_Os02g56610', 'LOC_Os01g61310', 'LOC_Os05g39500',
          'LOC_Os05g28040')
os_with_alog <- rbind(os, 
                      data.table(`Protein ID` = ALOG,
                                 Family = "ALOG"))

# combine families
combined_tfdb <- merge(os_with_alog,
                       families_raw[, .(
                           `Protein ID` = Gene_ID,
                           `planttfdb` = Family)],
                       all = TRUE)

# add class info for AP2s
ap2_family <- "AP2-EREBP"
ap2_classes <- c("ERF", "RAV", "AP2")
combined_tfdb[planttfdb %in% ap2_classes, Family := ap2_family]
combined_tfdb[planttfdb %in% ap2_classes, Class := planttfdb]

# add class info for HBs
hb_family <- "HB"
hb_classes <- c("TALE", "HD-ZIP", "WOX", "HB-other", "HB-PHD")
combined_tfdb[planttfdb %in% hb_classes, Family := hb_family]
combined_tfdb[planttfdb %in% hb_classes, Class := planttfdb]
combined_tfdb[Family == hb_family & is.na(planttfdb), Class := "HB-other"]

# add class info for MADS
mads_family <- "MADS"
mads_classes <- c("MIKC_MADS", "M-type_MADS")
combined_tfdb[planttfdb %in% mads_classes, Family := mads_family]
combined_tfdb[planttfdb %in% mads_classes, Class := planttfdb]

# match unnanotated Family to planttfdb
all_fam <- sort(combined_tfdb[, unique(Family)])
combined_tfdb[planttfdb == "MYB_related", planttfdb := "MYB-related"]
combined_tfdb[is.na(Family) & planttfdb %in% all_fam, Family := planttfdb]

# tidy up
combined_tfdb[, planttfdb := NULL]
tidy_tfdb <- combined_tfdb[!is.na(Family)]

# output
saveRDS(tidy_tfdb, output_file)

# log
sessionInfo()



