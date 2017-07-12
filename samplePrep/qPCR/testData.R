library(ggplot2)

meanprice <- tapply(diamonds$price, diamonds$cut, mean)
se <- tapply(diamonds$price, diamonds$cut, function(x) sqrt(var(x)/length(x))) # standard error of the mean
cut <- factor(levels(diamonds$cut), levels = levels(diamonds$cut))

g <- ggplot(data.frame(cut, meanprice, se),
            aes(x = cut, y = meanprice, ymin = meanprice - se, ymax = meanprice + se)) + 
  geom_bar(stat = 'identity', fill = 'grey', colour = 'black') +
  geom_errorbar(width = 0.2) +
  theme_minimal(base_size = 12)

g <- g + ylab((expression(atop(Normalised~relative~expression,(2^{- Delta * Delta * C[t]} %+-%~"95%"~CI*","~italic(n) == 3))))) + xlab(NULL)

cuts <- levels(diamonds$cut)
pvals <- sapply(1:length(cuts), function(i) 
  t.test(x = diamonds[diamonds$cut == cuts[1],'price'],
         y = diamonds[diamonds$cut == cuts[i],'price'],
         alternative = 'two.sided', conf.level = 0.95)$p.value)

labels <- data.frame(pvals, label = '', stringsAsFactors = FALSE)
labels[labels$pvals >= 0.05,'label'] <- 3
labels[labels$pvals < 0.05,'label'] <- 1
labels[labels$pvals < 0.005,'label'] <- 2
labels$asterisks <- c('*', '**', '')[as.numeric(labels$label)]

cairo_pdf('testOut.pdf', height = 7.5, width = 10, family = 'Verdana')
g + geom_point(aes(y = meanprice + se + 100, shape = labels$label, size = labels$label)) +
  scale_shape_manual(values = c('*',"#",NA),
                     labels = c(expression(paste(0.05 > italic(p), phantom() >= 0.005)), expression(italic(p) < 0.005), ''), guide = guide_legend(title=NULL, label.hjust = 0)) +
  scale_size_manual(values = c(8,5,0),
                    labels = c(expression(paste(0.05 > italic(p), phantom() >= 0.005)), expression(italic(p) < 0.005), ''), guide = guide_legend(title=NULL, label.hjust = 0)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1))
dev.off()

## not run
cairo_pdf('testOut2.pdf', height = 7.5, width = 10, family = 'Verdana')
g + geom_text(aes(y = meanprice + se + 50, label = labels$asterisks)) +
  geom_point(aes(shape = labels$label), colour = NA) +
  scale_shape_manual(values = c('1',"2",'3'),
                     labels = c(expression(paste('*'~~0.05 > italic(p), phantom() >= 0.005)), expression("**"~~italic(p) < 0.005), ''), guide = guide_legend(title=NULL, label.hjust = 0)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1))
dev.off()

g + geom_point(aes(y = meanprice + 150, shape = labels$label, size = labels$label)) +
  scale_shape_manual(values = c('*', '#', NA),
                     labels = c(expression(paste(0.05 > italic(p), phantom() >= 0.005)), expression(italic(p) < 0.005), ''), name = expression(atop("Student's"~italic(t)*"-test"~italic(vs)*". control,", italic(n) == 3))) +
  scale_size_manual(values = c(8,5,0),
                    labels = c(expression(paste(0.05 > italic(p), phantom() >= 0.005)), expression(italic(p) < 0.005), ''), name = expression(atop("Student's"~italic(t)*"-test"~italic(vs)*". control,", italic(n) == 3)))
