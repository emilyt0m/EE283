library( "DESeq2" )
sampleInfo = read.table("shortRNAseq.txt", header=TRUE)
sampleInfo$FullSampleName = as.character(sampleInfo$FullSampleName)

countdata = read.table("fly_counts2.txt", header=TRUE, row.names=1)
countdata = countdata[ ,6:ncol(countdata)]

temp = colnames(countdata)
temp = gsub("RNAseq.bam.","",temp)
temp = gsub(".bam","",temp)
colnames(countdata) = temp

cbind(temp,sampleInfo$FullSampleName,temp == sampleInfo$FullSampleName)

dds = DESeqDataSetFromMatrix(countData=countdata, colData=sampleInfo, design=~TissueCode)
dds <- DESeq(dds)
res <- results( dds )
res

plotMA( res, ylim = c(-1, 1) )
plotDispEsts( dds )
hist( res$pvalue, breaks=20, col="grey" )

rld = rlog( dds )
head( assay(rld) )
mydata = assay(rld)

sampleDists = dist( t( assay(rld) ) )
sampleDistMatrix = as.matrix( sampleDists )
rownames(sampleDistMatrix) = rld$TissueCode
colnames(sampleDistMatrix) = NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

print( plotPCA( rld, intgroup = "TissueCode") )

library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", margins=c(6,6), col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

library(EnhancedVolcano)
res <- res[order(res$padj),]
EnhancedVolcano(res, lab= rownames(res), x = 'log2FoldChange', y = 'pvalue', pointSize = 1.0, labSize = 2.0, legendLabels=c('Not sig.','Log (base 2) FC','p-value','p-value & Log (base 2) FC'), legendPosition = 'right', legendLabSize = 8, legendIconSize = 2.0)
