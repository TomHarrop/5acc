#!/usr/bin/env Rscript

library(data.table)

#############################
# GRanges wrapper functions #
#############################

FormatRegionsAsGRange <- function(range.data.table, seqinfo = seq.info){
  # start with unique start/end coordinates
  unique.ranges <- unique(range.data.table[, .(
    CHR, `window-start`, `window-end`)])
  # make GRanges object
  my.gr <- GenomicRanges::makeGRangesFromDataFrame(
    df = unique.ranges,
    start.field = "window-start",
    end.field = "window-end",
    seqnames.field = "CHR",
    ignore.strand = TRUE,
    keep.extra.columns = FALSE,
    seqinfo = seqinfo)
  # merge overlapping or adjacent windows
  GenomicRanges::reduce(my.gr)
}

AddTestMetadataColumn <- function(test.name, granges.list){
  # add test.name as metadata column
  my.gr <- granges.list[[test.name]]
  GenomicRanges::mcols(my.gr)$test <- test.name
  my.gr
}

ExtractGenesInTestWindows <- function(test.range, gtf.genes = gtf.genes){
  # find overlaps
  query.results <- as.data.table(
    GenomicRanges::findOverlaps(query = test.range,
                                subject = gtf.genes))
  # add region info
  query.results[, c("seqnames", "start", "end", "width", "strand", "test") := 
                  as.data.frame(test.range[queryHits]), by = queryHits]
  
  # add gene names
  query.results[, gene := names(gtf.genes)[subjectHits], by = subjectHits]
  
  # return
  query.results[, .(
    gene, Chr = seqnames, window.start = start, window.end = end, test
  )]
}

########################
# Get seqinfo from GTF #
########################

gtf.file <- 
  "data/genome/os/Osativa_323_v7.0.gene_exons.cuffcomp.rRNAremoved.gtf"

gtf <- rtracklayer::import.gff(gtf.file,
                               format = "gtf",
                               genome = "Osativa_323_v7.0",
                               feature.type = "exon")

seq.info <- GenomicRanges::seqinfo(gtf)

# split gtf by gene
gtf.split <- GenomicRanges::split(
  gtf, GenomicRanges::elementMetadata(gtf)$gene_name)
gtf.genes <- GenomicRanges::reduce(gtf.split)

#########################################
# Process test results into GRangesList #
#########################################

test.result.files <- list.files(
  "data/genes_under_selection",
  pattern = ".csv",
  full.names = TRUE)

names(test.result.files) <- gsub(".+_[[:blank:]]?([[:alnum:]]+).+",
                                 "\\1",
                                 x = basename(test.result.files))

test.results <- lapply(test.result.files,
                       FUN = fread,
                       skip = "CHR",
                       encoding = "UTF-8")

test.ranges <- lapply(test.results,
                      FUN = FormatRegionsAsGRange,
                      seqinfo = seq.info)

test.ranges.with.metadata <- lapply(names(test.ranges),
                                    FUN = AddTestMetadataColumn,
                                    granges.list = test.ranges)
names(test.ranges.with.metadata) <- names(test.ranges)

#######################################
# find genes overlapping test regions #
#######################################

window.genes.list <- lapply(test.ranges.with.metadata,
                            FUN = ExtractGenesInTestWindows,
                            gtf.genes = gtf.genes)
window.genes <- rbindlist(window.genes.list)
setkey(window.genes, gene)

###################
# add annotations #
###################

annotations <- window.genes[, oryzr::LocToGeneName(gene)]
setkey(annotations, MsuID)
annotated.window.genes <- annotations[window.genes, .(
  Chr, window.start, window.end, test,
  gene = MsuID, symbol = symbols, annotation = MsuAnnotation,
  OgroObjective, OgroRef
)]

###############
# Save output #
###############

outdir = "output/genes_under_selection"
if(!dir.exists(outdir)){
  dir.create(outdir)
}
saveRDS(annotated.window.genes, paste0(outdir, "/annotated_window_genes.RDS"))
write.table(annotated.window.genes,
            paste0(outdir, "/annotated_window_genes.tab"),
            sep = "\t",
            quote = FALSE,
            na = "",
            row.names = FALSE)

sinfo <- rutils::GitSessionInfo()
writeLines(sinfo, paste0(outdir, "/SessionInfo.genes_under_selection.txt")
