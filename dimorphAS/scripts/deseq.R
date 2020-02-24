### DESeq ###
library(DESeq2)
library(RColorBrewer)
library(pheatmap)

# input varibales
tissue<-'Brain-Cortex' # variable
mcount<-76
fcount<-37


setwd("/projects/kesara/Collaborator/Peter")
ifile<-paste(c('/projects/kesara/Collaborator/Peter/',tissue,'.rbp_genes_count.csv'),  sep="", collapse='')
count <- read.csv(ifile, row.names= 1, header = T)
head(count)

countdata <- as.matrix(count[,2:ncol(count)])
head(countdata)
dim(countdata)

# Assign conditions
sex <- factor(c(rep("M", mcount), rep("F", fcount)))
#sex

# DESeq
# Create a coldata frame and instantiate the DESeqDataSet. See ins?DESeqDataSetFromMatrix
coldata<-data.frame(row.names=colnames(countdata), sex)
#coldata
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design =~ sex) # RK
coldata

# Pre-filtering; sum of read count
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1, ]
nrow(dds)

# rlog
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

# PCA
ofile<-paste(c(tissue,'.RBP.PCA.png'), sep="", collapse='')
png(ofile, 1000, 1000, pointsize=20)
plotPCA(rld, intgroup = "sex")
dev.off()

# sample distance
sampleDists <- dist(t(assay(rld)))
#sampleDists

# Heatmap of sample-to-sample distances using the rlog-transformed values.
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(coldata$sex, colnames(countdata), sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

ofile<-paste(c(tissue,'.RBP.pheatmap.png'), sep="", collapse='')
png(ofile, 1000, 1000, pointsize=20)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists, col = colors)
dev.off()

# Run the DESeq pipeline; by default test is 'Wald test'
ofile<-paste(c(tissue,'.RBP.diff_exp.csv'),  sep="", collapse='')
dds <- DESeq(dds)
res <- results(dds)
sig_res<-res$padj<0.05
sig_res
#table(res$padj<0.05)
#res <- res[order(res$padj), ]
write.csv(res, file=ofile)