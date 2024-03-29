#Load counts table
c_vs_d_countsTable <- read.csv("c vs. d counts.csv", header=T, stringsAsFactors=FALSE)

rownames(c_vs_d_countsTable) <- c_vs_d_countsTable[,1]
c_vs_d_countsTable <- c_vs_d_countsTable[,-1]

library(DESeq2)
c_vs_d_conds <- factor(c("CTRL", "CTRL", "KO_D","KO_D"))
c_vs_d_colData = data.frame(condition = c_vs_d_conds)
dds <- DESeqDataSetFromMatrix(c_vs_d_countsTable, c_vs_d_colData, design=~ condition)

# Exploratory data analysis
# Relationship between variance (sigma^2) and mean (mu) of genes across samples
mu = rowMeans(c_vs_d_countsTable)
sigma2 = apply(c_vs_d_countsTable,1,var)
plot(log(mu), log(sigma2), xlim = c(0, log(max(mu))), ylim = c(0, log(max(mu))), pch=16, cex=0.3)

abline(0, 1, col="red")

#Correct for heteroskedasticity
rld <- rlog(dds)
colnames(rld) = colnames(c_vs_d_countsTable)

plotPCA(rld, intgroup = "condition")

vsd <- varianceStabilizingTransformation(dds)
colnames(vsd) = colnames(knockoutcountsTable)

plotPCA(vsd, intgroup = "condition")

#Differential expression analysis
c_vs_d_dds <- DESeq(dds)
c_vs_d_res <- results(c_vs_d_dds)
#Meaning of columns
mcols(c_vs_d_res, use.names=TRUE)

summary(c_vs_d_res)

#Sort results by adjusted p value in ascending order
c_vs_d_res <- c_vs_d_res[order(res$padj),]

#Subset results for genes with BH adjusted p value < 0.05
c_vs_d_deg = subset(c_vs_d_res, padj < 0.05)

write.csv(as.data.frame(c_vs_d_res), file = "c vs. d degresults.csv")

#Gene clusters and heatmap
library("RColorBrewer")
library("gplots")
ramp <- 1:3/3
cols <- c(rgb(ramp,0,0), rgb(0,ramp,0), rgb(0,0,ramp), rgb(ramp,0,ramp))
topVarGenes <- head(order(genefilter::rowVars(assay(rld)), decreasing=T), 25)
heatmap.2(assay(rld)[topVarGenes,], scale="row", trace="none", density.info="none", dendrogram="none", margins=c(6,6), colv="FALSE", col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

#Volcano Plot
library(EnhancedVolcano)
EnhancedVolcano(c_vs_d_res, lab= rownames(c_vs_d_res), x = 'log2FoldChange', y = 'pvalue')
