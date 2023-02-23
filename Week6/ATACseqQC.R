install.packages("tidyverse")
install.packages("BiocManager")
BiocManager::install("ATACseqQC")

library(ATACseqQC)

bamfile <- c("A4_ED_2.dedup.bam", "A4_ED_3.dedup.bam", "A4_ED_4.dedup.bam", "A4_WD_1.dedup.bam","A4_WD_2.dedup.bam", "A4_WD_4.dedup.bam")
bamfile.labels <- gsub(".bam", "", basename(bamfile))
fragSize <- fragSizeDist(bamfile, bamfile.labels)

estimateLibComplexity(readsDupFreq(bamfile[1]))
estimateLibComplexity(readsDupFreq(bamfile[2]))
estimateLibComplexity(readsDupFreq(bamfile[3]))
estimateLibComplexity(readsDupFreq(bamfile[4]))
estimateLibComplexity(readsDupFreq(bamfile[5]))
estimateLibComplexity(readsDupFreq(bamfile[6]))

outPath <- "A4_WD_4 splitBam"
dir.create(outPath)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
library(BSgenome.Dmelanogaster.UCSC.dm6)

possibleTag <- combn(LETTERS, 2)
possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]),
                 paste0(possibleTag[2, ], possibleTag[1, ]))
library(Rsamtools)
bamTop100 <- scanBam((bamfile[6]), yieldSize = 100, param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]

gal <- readBamFile(bamfile[6], tag=tags, asMates=TRUE)
gal1 <- shiftGAlignmentsList(gal)
shiftedBamfile <- file.path(outPath, "A4_WD_4 shifted.bam")
export(gal1, shiftedBamfile)
                      