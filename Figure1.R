### Figure 1
# load the following objects from https://osf.io/c84vs/
dds_mrg <- readRDS("dds_mrg.Rds")


### Figure 1C
# PCA plot of normal diff for Ob and Ad seperately
library(DESeq2)
library(genefilter)
library(ggplot2)

# Get rlog values of samples of interest and subject to PCA
rld_RNA <- data.frame(assay(vst(dds_mrg, blind=F)))
rld_RNA <- rld_RNA[,c(grep("_1",colnames(rld_RNA)),grep("_2",colnames(rld_RNA)))]
colnames(rld_RNA)
rld_Ob <- rld_RNA[,c(1:6,18:23)]
rld_Ad <- rld_RNA[,c(6:11,23:28)]

# PCA plot
# Osteoblasts
object <- rld_Ob
intgroup <- c('Group','Replicate')
intgroup.df <- data.frame(colData(dds_mrg))
intgroup.df <- intgroup.df[colnames(object),]
ntop=500

rv <- rowVars(object)
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(object[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
if (!all(intgroup %in% names(intgroup.df))) {
  stop("the argument 'intgroup' should specify columns of intgroup.df")
}
rownames(intgroup.df) <- intgroup.df$Sample
group <- if (length(intgroup) > 1) {
  factor(apply( intgroup.df, 1, paste, collapse=" : "))
} else {
  intgroup.df[[intgroup]]
}
dat <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], PC4=pca$x[,4], group=group, intgroup.df, name=colnames(object))
attr(dat, "percentVar") <- percentVar[1:4]

ratio.values <- (max(dat$PC1)-min(dat$PC1))/(max(dat$PC2)-min(dat$PC2))
ggplot(data=dat, aes_string(x="PC1", y="PC2", color="Group"))+ geom_point(size=3) + coord_fixed(ratio.values) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance"))

# Adipocytes
object <- rld_Ad
intgroup <- c('Group','Replicate')
intgroup.df <- data.frame(colData(dds_mrg))
intgroup.df <- intgroup.df[colnames(object),]
ntop=500

rv <- rowVars(object)
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(object[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
if (!all(intgroup %in% names(intgroup.df))) {
  stop("the argument 'intgroup' should specify columns of intgroup.df")
}
rownames(intgroup.df) <- intgroup.df$Sample
group <- if (length(intgroup) > 1) {
  factor(apply( intgroup.df, 1, paste, collapse=" : "))
} else {
  intgroup.df[[intgroup]]
}
dat <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], PC4=pca$x[,4], group=group, intgroup.df, name=colnames(object))
attr(dat, "percentVar") <- percentVar[1:4]

ratio.values <- (max(dat$PC1)-min(dat$PC1))/(max(dat$PC2)-min(dat$PC2))
ggplot(data=dat, aes_string(x="PC1", y="PC2", color="Group"))+ geom_point(size=3) + coord_fixed(ratio.values) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance"))

# Final aesthetics, such as colors and line drawing, were done in Illustrator