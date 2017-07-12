library(ggplot2)
library(plyr)
library(reshape2)

setwd('/home/tom/Desktop/5accessions/samplePrep/testLibs/markdup/')

picardOutput <- dir(pattern = '.ISM')

readHistData <- function(picardFile){
  return(read.table(picardFile, quote = '#', skip = 11, header = TRUE, stringsAsFactors = FALSE))}

histData <- lapply(picardOutput, readHistData)
names(histData) <- unlist(strsplit(picardOutput, '\\.'))[c(1,3)]
for (i in 1:length(histData)){histData[[i]]$name <- names(histData)[i]}

histData.frame <- do.call(rbind, histData)

plotData <- melt(histData.frame, id.vars = c('name', 'insert_size'), measure.vars = c('All_Reads.fr_count', 'All_Reads.rf_count'))

g <- ggplot(plotData, aes(x = insert_size, y = value, colour = variable)) +
#  geom_line() +
  theme_minimal() +
  scale_colour_brewer(palette = 'Set1') + 
  stat_smooth(method = 'loess', span = 0.05, se = FALSE) +
  facet_grid(.~name)
#  geom_vline(x = c(145-59,145,145+59))

gb <- ggplot_build(g)$data[[1]]
#gb$sample <- as.factor(c('RF', 'FR', 'FR', 'RF')[gb$group])
gb$sample <- as.factor(c('FR', 'FR', 'RF', 'RF')[gb$group])

ggplot(gb, aes(x=x, y=y, fill = sample, colour = sample)) + geom_area(alpha = 0.6) +
  scale_fill_brewer(palette = 'Set1') + geom_line(colour = 'black')

as.factor(c('RF', 'FR', 'FR', 'RF')[gb$group])
as.factor(c('RF', 'FR', 'RF', 'FR')[gb$group])

###############
### NOT RUN ###
###############

ggplot(histData.frame, aes(x = insert_size, y = All_Reads.fr_count)) + stat_smooth()

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
