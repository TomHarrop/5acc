library(ggplot2)
library(plyr)

setwd('/home/tom/Desktop/5accessions/samplePrep/testLibs/bwamem')

SAMs <- dir(pattern = '.sam')

sizeMetrics <- function(SAMfile) {
  
  # Calculate the number of columns in the body of the SAM file
  
  SAMbody <- read.table(SAMfile, header = FALSE, comment.char = '@', sep = '\t', stringsAsFactors = FALSE, fill = TRUE, nrows = 5000)
  SAMwidth <- dim(SAMbody)[2]
  
  # Set up col.names for SAM file import
  
  SAMnames<- c('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'MRNM', 'MPOS', 'ISIZE', 'SEQ', 'QUAL', 'TAG', 'VTYPE', 'VALUE')
  
  if (length(SAMnames) == SAMwidth) {
    SAM_colnames <- SAMnames
  } else if (length(SAMnames) < SAMwidth) {
    i <- SAMwidth-length(SAMnames)
    newCols <- seq(SAMwidth + 1, SAMwidth + i)
    newCols <- sapply(newCols, function(x) paste("V", x, sep = ''))
    SAM_colnames <- c(SAMnames, newCols)
  } else if (length(SAMnames) > SAMwidth) {
    SAM_colnames <- SAMnames[1:SAMwidth]
  }
  
  # Add the colnames to the SAMFile
  
  colnames(SAMbody) <- SAM_colnames
  SAMbody$ISIZE <- as.numeric(SAMbody$ISIZE)
  
  # Calculate the mean insert size from high-quality mappings. Exclude inserts >
  # 500 bp
  
  matePairs <- subset(SAMbody, MRNM == '=' & MAPQ == '60')
  matePairsPosLength <- subset(matePairs, ISIZE > 0 & ISIZE < 500)
  matePairsPosLength$samfile <- SAMfile
  return(matePairsPosLength)
}

SAMmetrics <- lapply(SAMs, function(x) sizeMetrics(x)[,c('samfile',  'ISIZE', 'QNAME')])

SAMmetrics.frame <- do.call(rbind, SAMmetrics)
SAMmetrics.frame$samfile <- revalue(factor(SAMmetrics.frame$samfile), c("ob.bwamem.sam"="O. barthii\nN1R2", "osn.bwamem.sam"="O. sativa Nipponbare\nN3R3"))
cdf <- ddply(SAMmetrics.frame, "samfile", summarise, ISIZE.mean=mean(ISIZE), ISIZE.sd = sd(ISIZE))
cdf$label <- paste('mean = ', round(cdf$ISIZE.mean, 0), '\nSD = ', round(cdf$ISIZE.sd, 0), sep = '')

g <- ggplot(SAMmetrics.frame, aes(x =ISIZE)) + theme_minimal(base_size = 12) +
  geom_density(aes(fill = samfile), colour = NA) + facet_wrap(~ samfile) +
  scale_fill_brewer(palette = "Set1", guide = FALSE) +
  geom_line(stat = 'density', size = 0.5) +
  geom_text(data = cdf, aes(x = 250, y = 0.0075, label = label), size = 4) +
  xlab(NULL) + ylab(NULL) + ggtitle("Insert Size")
cairo_pdf('insertSize.pdf', family = 'Inconsolata', width = 10, height = 7.5)
g
dev.off()
