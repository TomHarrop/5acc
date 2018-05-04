#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(valr)

###########
# GLOBALS #
###########

os_gff_file <- snakemake@input[["os_gff_file"]]
os_gtf_file <- snakemake@input[["os_gtf_file"]]
seqlengths_file <- snakemake@input[["seqlengths_file"]]
irgsp_gff_file <- snakemake@input[["irgsp_gff_file"]]
osa1r7_gff_file <- snakemake@input[["osa1r7_gff_file"]]
osa1_mirbase_gff_file <- snakemake@input[["osa1_mirbase_gff_file"]]
tigr_repeats_fa <- snakemake@input[["tigr_repeats_fa"]]

star_index_dir <- snakemake@params[["star_index_dir"]]
cpus <- snakemake@threads[[1]]
log_file <- snakemake@log[["log"]]

shuffled_gtf_file <- snakemake@output[["shuffled_gtf"]]

# dev
# cpus <- 8
# os_gff_file <- "data/genome/os/Osativa_323_v7.0.gene_exons.gff3"
# os_gtf_file <- "output/010_data/Osativa_323_v7.0.gene_exons.cuffcomp.rRNAremoved.gtf"
# star_index_dir <- "output/010_data/star-index"
# seqlengths_file <- "output/010_data/star-index/chrNameLength.txt"
# irgsp_gff_file <- "data/genome/os/irgsp1_rRNA_tRNA.gff"
# osa1r7_gff_file <- "data/genome/os/rice_osa1r7_rm.gff3"
# osa1_mirbase_gff_file <- "data/genome/os/osa.gff3"
# tigr_repeats_fa <- "data/genome/os/TIGR_Oryza_Repeats.v3.3_0_0.fsa"

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# load tbl genome
genome <- read_genome(seqlengths_file)

# get genes
os_gff_genes <- import.gff(os_gff_file, feature.type = "gene")

# slop genes
slopped_genes_tbl <- bed_slop(as.tbl_interval(os_gff_genes),
                              genome,
                              both = 100)

# fix irgsp gff
irgsp_tmp1 <- tempfile(fileext = ".gff")
irgsp_tmp2 <- tempfile(fileext = ".gff")
irgsp_tmp3 <- tempfile(fileext = ".gff")
system2("sed",
        args = c("282d", irgsp_gff_file),
        stdout = irgsp_tmp1,
        stderr = log_file)
system2("sed",
        args = c("537d", irgsp_tmp1),
        stdout = irgsp_tmp2,
        stderr = log_file)
system2("sed",
        args = c("913d", irgsp_tmp2),
        stdout = irgsp_tmp3,
        stderr = log_file)
irgsp_gff <- import.gff(irgsp_tmp3)

# rename chromosomesq
slr <- gsub("0(\\d)",
            "\\1",
            sub("chr",
                "Chr",
                seqlevels(irgsp_gff)))
names(slr) <- seqlevels(irgsp_gff)
seqlevels(irgsp_gff) <- slr
irgsp_tbl <- as.tbl_interval(irgsp_gff)

# load osa1r7 gff
osa1r7_gff <- import.gff(osa1r7_gff_file)
osa1r7_tbl <- as.tbl_interval(osa1r7_gff)

# load osa.gff3 miRBase miRNAs
osa1_mirbase_gff <- import.gff(osa1_mirbase_gff_file)
osa1_mirbase_tbl <- as.tbl_interval(osa1_mirbase_gff)

# simulate reads from tigr repeats
wgsim1 <- tempfile(fileext = ".wgsim.1.fq")
wgsim2 <- tempfile(fileext = ".wgsim.2.fq")
system2("wgsim",
        args = c("-e", "0",
                 "-1", "55",
                 "-2", "55",
                 "-r", "0",
                 "-R", "0",
                 "-X", "0",
                 "-d", "0",
                 "-s", "0",
                 tigr_repeats_fa,
                 wgsim1,
                 wgsim2),,
        stdout = log_file
        stderr = log_file)

# map tigr repeats
star_outdir <- tempdir()
prefix <- paste(star_outdir, "TIGR_Oryza_Repeats.", sep = "/")
system2("STAR",
        args = c("--runThreadN",
                 cpus,
                 "--genomeDir",
                 star_index_dir,
                 "--outSAMtype",
                 "BAM SortedByCoordinate",
                 "--outFilterMultimapNmax",
                 "-1",
                 "--outBAMcompression",
                 "10 ",
                 "--readFilesIn",
                 wgsim1, wgsim2,
                 "--outFileNamePrefix",
                 prefix),,
        stdout = log_file
        stderr = log_file)

# convert BAM to bed
rpt_bed <- tempfile(fileext = ".bed6")
star_bamfile <- paste0(prefix, "Aligned.sortedByCoord.out.bam")
system2("bedtools",
        args = c("bamtobed",
                 "-i", star_bamfile),
        stdout = rpt_bed,
        stderr = log_file)

# read bed hits
rpt_hits <- import.bed(rpt_bed)
rpt_tbl <- as.tbl_interval(rpt_hits)

# merge with valr
all_tbl <- dplyr::bind_rows(slopped_genes_tbl,
                            irgsp_tbl,
                            osa1r7_tbl,
                            osa1_mirbase_tbl,
                            rpt_tbl)
merged_tbl <- bed_merge(all_tbl)

# prepare a dummy GFF for shuffling
os_gff_exons <- import.gff(os_gtf_file, feature.type = "exon", format = "gtf")
grl <- GenomicRanges::reduce(split(os_gff_exons,
                                   elementMetadata(os_gff_exons)$gene_name))
gtf_reduced <- unlist(os_gff_exons, use.names = FALSE)

# add metadata
elementMetadata(gtf_reduced)$widths <- width(gtf_reduced)

# calculate feature lengths with dplyr
feature_length_tbl <- group_by(as.data.frame(gtf_reduced), gene_name) %>%
  summarize(length = sum(widths))
feature_lengths <- data.table(Length = feature_length_tbl$length,
                              rn = feature_length_tbl$gene_name,
                              key = "rn")
to_shuffle <- feature_lengths[Length < quantile(feature_lengths$Length, 0.9),
                              unique(rn)]

# generate dummy ranges
gene_chromosome <- unique(
  data.table(rn = os_gff_genes$Name,
             seqid = as.character(GenomeInfoDb::seqnames(os_gff_genes)),
             strand = as.character(rtracklayer::strand(os_gff_genes)),
             key = "rn"))
dummy_gff_dt <- gene_chromosome[feature_lengths, .(
  chrom = seqid,
  source = 'phytozomev10',
  type = 'CDS',
  start = 1,
  end = Length,
  score = ".",
  strand,
  phase = ".",
  ID = rn
)]
dummy_gff_tbl <- as.tbl_interval(dummy_gff_dt)

# shuffle
shuffled_gtf <- bed_shuffle(
  dummy_gff_tbl %>% filter(
    (!chrom %in% c("ChrSy", "ChrUn")) &
      ID %in% to_shuffle),
  genome,
  excl = merged_tbl,
  within = TRUE,
  seed = 1)

# convert to Granges
shuffled_gr <- makeGRangesFromDataFrame(shuffled_gtf,
                                        keep.extra.columns=FALSE,
                                        ignore.strand=TRUE,
                                        seqinfo=NULL,
                                        seqnames.field="chrom",
                                        start.field="start",
                                        end.field="end",
                                        strand.field="strand")
names(shuffled_gr) <- shuffled_gtf$ID

# write output
export(shuffled_gr, shuffled_gtf_file, "gff3")

# write log
sessionInfo()
