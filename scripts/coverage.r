args<-commandArgs(TRUE)

all_results <- read.table(args[1], header=TRUE, row.names=1)
genes <- read.table(args[2])
results <- all_results[,grep(sprintf("above_%s$", args[3]),names(all_results))]
mat <- as.matrix(results)
results$means <- rowMeans(mat)
results$gene <- genes[,1]
results$coords = rownames(results)

library(ggplot2)
library(plyr)

pdf(args[4])
d_ply(results, .var = "gene", function(x) {
p <- ggplot( x, aes(y=means, x=factor(coords))) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(ylab(sprintf("percentage exon at least %sX", args[3]))) + ggtitle(unique(x$gene)) + ylim(0, 100)
print(p)
})
dev.off()
