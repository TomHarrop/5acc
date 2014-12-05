setwd('/home/tom/Desktop/5accessions/samplePrep/qPCR')

nanoQuantData <- readRDS('../nanoQuant/nanoQuandData.RDS')
nanoQuantData$NIP.concentration <- c(79.0,376.8,389.6,357.9,347.0,309.8,586.6,411.5,709.0,457.1,530.1,510.9)
#update re-runs
nanoQuantData['N1R1','NIP.concentration'] <- 148.0
nanoQuantData['N3R1','IR64.concentration'] <- 94.7

stageFactor <- factor(rep(c('N1', 'N2', 'N3', 'N4'), each = 3))

maxTable <- apply(nanoQuantData, 2, function(x) by(x, stageFactor, which.max))

columnName <- 'IR64.concentration'

getMaxValues <- function(columnName, nanoQuantData, maxTable){
  indexes <- maxTable[,columnName] + c(0, 3, 6, 9)
  retValues <- c(nanoQuantData[indexes,columnName])
  names(retValues) <- rownames(nanoQuantData)[indexes]
  return(retValues)
}

getMaxValues('NIP.concentration', nanoQuantData, maxTable)

species <- factor(rep(unique(unlist(lapply(colnames(nanoQuantData), function(x) strsplit(x, split = '\\.')[[1]][1]))), each = 4), levels = c('NIP', 'IR64', 'OR', 'OG', 'OB'))
cdnaTable <- data.frame(species = species[order(species)], stage = factor(rep(c('N1', 'N2', 'N3', 'N4'), 5)))
cdnaTable$sample <- 0
cdnaTable$concentration <- 0
cdnaTable$volume <- 0
cdnaTable$water <- 0
cdnaTable$amount <- 0
cdnaTable$totalamount <- 0

for (x in levels(species)){  
  cdnaTable[cdnaTable$species == x,'concentration'] <- getMaxValues(paste(x, '.concentration', sep = ''), nanoQuantData, maxTable)
  cdnaTable[cdnaTable$species == x,'sample'] <- names(getMaxValues(paste(x, '.concentration', sep = ''), nanoQuantData, maxTable))
  cdnaTable[cdnaTable$species == x,'totalamount'] <- min(getMaxValues(paste(x, '.concentration', sep = ''), nanoQuantData, maxTable)) * 11
  cdnaTable[cdnaTable$species == x,'amount'] <- round(cdnaTable[cdnaTable$species == x,'totalamount'] / 4, 1)
  cdnaTable[cdnaTable$species == x,'volume'] <- round(cdnaTable[cdnaTable$species == x,'amount'] / cdnaTable[cdnaTable$species == x,'concentration'], 1)
  cdnaTable[cdnaTable$species == x,'totalamount'] <- round(cdnaTable[cdnaTable$species == x,'totalamount']/1000, 1)
  cdnaTable[cdnaTable$species == x,'water'] <- 11-sum(cdnaTable[cdnaTable$species == x,'volume'])
}

write.csv(cdnaTable, file = 'cdnaTable.csv', row.names = FALSE)
