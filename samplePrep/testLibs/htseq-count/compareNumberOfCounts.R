setwd('/home/tom/Desktop/5accessions/samplePrep/testLibs/htseq-count')

library(DESeq2)

fileName <- dir('.', pattern = 'htseq-count$', include.dirs = FALSE)
sampleName <- sapply(fileName, function(x) toupper(strsplit(x, split = '\\.')[[1]][1]))
sampleTable <- data.frame(sampleName, fileName)
rownames(sampleTable) <- NULL

parameter <- factor(sampleTable$sampleName)
sampleTable$parameter <- parameter
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = '.', design = ~ parameter)

colSums(counts(dds))
counts <- counts(dds)

which.max(counts[,'SSREV'])
counts['LOC_Os09g08910',]

which.max(abs(counts[,'FR-FIRSTSTRAND-SWAPPED'] - counts[,'FR-SECONDSTRAND-SWAPPED']))
counts['LOC_Os02g38920',]

ordering <- abs(counts[,'FR-FIRSTSTRAND-SWAPPED'] - counts[,'FR-SECONDSTRAND-SWAPPED'])
head(counts[rev(order(ordering)),])

counts(dds)['LOC_Os03g11614',]
counts(dds)['LOC_Os10g33780',]

###############
### NOTE!!! ###
###############

### It appears the reads are incorrectly split by the fastq-dump command, so we 
### need to supply them to tophat in reverse order and use --library-type 
### fr-secondstrand, e.g. tophat [options] reads.2.fastq reads.1.fastq. The use
### of fr-secondstrand is based on visual inspection of the mappings in
### /home/tom/Data/zhang_data/subset/tophat, e.g. at LOC_Os02g38920 and
### LOC_Os12g44350. For how these regions were found see
### /home/tom/Data/zhang_data/subset/htseq-count/compareNumberOfCounts.R

