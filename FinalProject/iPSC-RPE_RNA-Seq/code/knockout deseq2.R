#Load counts table
knockoutcountsTable <- read.csv("knockoutcounts.csv", header=T, stringsAsFactors=FALSE)

rownames(knockoutcountsTable) <- knockoutcountsTable[,1]
knockoutcountsTable <- knockoutcountsTable[,-1]

library(DESeq2)
knockoutconds <- factor(c("KO1","KO1","KO2","KO2"))
knockoutcolData = data.frame(condition = knockoutconds)
dds <- DESeqDataSetFromMatrix(knockoutcountsTable, knockoutcolData, design=~ condition)

# Exploratory data analysis
# Relationship between variance (sigma^2) and mean (mu) of genes across samples
mu = rowMeans(knockoutcountsTable)
sigma2 = apply(knockoutcountsTable,1,var)
plot(log(mu), log(sigma2), xlim = c(0, log(max(mu))), ylim = c(0, log(max(mu))), pch=16, cex=0.3)

abline(0, 1, col="red")

#Correct for heteroskedasticity
rld <- rlog(dds)
colnames(rld) = colnames(knockoutcountsTable)

plotPCA(rld, intgroup = "condition")

vsd <- varianceStabilizingTransformation(dds)
colnames(vsd) = colnames(knockoutcountsTable)

plotPCA(vsd, intgroup = "condition")

#Differential expression analysis
knockoutdds <- DESeq(dds)
knockoutres <- results(knockoutdds)
#Meaning of columns
mcols(knockoutres, use.names=TRUE)

summary(knockoutres)

#Sort results by adjusted p value in ascending order
knockoutres <- knockoutres[order(res$padj),]

#Subset results for genes with BH adjusted p value < 0.05
knockoutdeg = subset(knockoutres, padj < 0.05)

write.csv(as.data.frame(knockoutres), file = "knockoutdegresults.csv")

#Gene clusters and heatmap
library("RColorBrewer")
library("gplots")
ramp <- 1:3/3
cols <- c(rgb(ramp,0,0), rgb(0,ramp,0), rgb(0,0,ramp), rgb(ramp,0,ramp))
topVarGenes <- head(order(genefilter::rowVars(assay(rld)), decreasing=T), 25)
heatmap.2(assay(rld)[topVarGenes,], scale="row", trace="none", density.info="none", dendrogram="none", margins=c(6,6), colv="FALSE", col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

#Volcano Plot
library(EnhancedVolcano)
EnhancedVolcano(knockoutres, lab= rownames(knockoutres), x = 'log2FoldChange', y = 'pvalue')
