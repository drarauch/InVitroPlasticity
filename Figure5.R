### Figure 5
# load the following objects (see XX for more information on how these were generated)
dds_mrg <- readRDS("dds_mrg.Rds")
RNA_tbl <- readRDS("RNA_tbl.Rds")

### Figure 5B
# PCA plot of normal diff up to day 7 and changing media at day 3 of Diff
library(DESeq2)
library(genefilter)
library(ggplot2)

# Get rlog values of samples of interest and subject to PCA
rld_RNA <- data.frame(assay(vst(dds_mrg, blind=F)))
rld_RNA <- rld_RNA[,c(grep("_1",colnames(rld_RNA)),grep("_2",colnames(rld_RNA)))]
colnames(rld_RNA)
rld_RNA <- rld_RNA[,c(2:10,12:15,19:27,29:32)]

# PCA plot
object <- rld_RNA
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

# Final aesthetics, such as color coding and line drawing, were done in Illustrator


### Figure 5C
# Scaled RNA-expression of normal diff up to day 3 and changing media at day 3 of Diff for cluster specific gene expression from Rauch et al., https://doi.org/10.1038/s41588-019-0359-1 
library(DESeq2)

# Get rlog values of samples of interest
rld_RNA <- data.frame(assay(vst(dds_mrg, blind=F)))
rld_RNA <- rld_RNA[,c(grep("_1",colnames(rld_RNA)),grep("_2",colnames(rld_RNA)))]
colnames(rld_RNA)
# order samples
rld_RNA <- rld_RNA[,c(14,12,3:9,13,15,31,29,20:26,30,32)]
# get average
rld_RNA <- (rld_RNA[,1:11] + rld_RNA[,12:22])/2
# scale data
rld_RNA <- data.frame(t(scale(t(rld_RNA))))

# boxplot the scaled signal
par(mfrow=c(1,3),pty='s')
for(i in c(5,9,1)) {
  genes <- TiCoNE[TiCoNE$tz_cluster_new == i,"RefSeqID"]
  genes <- genes[genes %in% rownames(rld_RNA)]
  tmp <- rld_RNA[genes,]
  Mycol <- rainbow(14)[i]
  MycolB <- 'blue'
  MycolR <- 'red'
  boxplot(tmp[,1],
          tmp[,2],
          tmp[,3],
          tmp[,4],
          tmp[,5],
          tmp[,6],
          tmp[,7],
          tmp[,8],
          tmp[,9],
          tmp[,10],
          tmp[,11],
          outline=F,col=c('blue','red','white','white','white','white','white','white','white','blue','red'),xlab="", ylab="", yaxt="n", xaxt="n")
  axis(2, at=c(-2,0,2), lab=c(-2,0,2), cex.axis=1)
  axis(1, at=c(1:9), lab=FALSE)
}


### Figure 5D
# Barplots for gene of interest (GOI)
library(ggpubr)
library(DESeq2)

# Define GOI
GOI <- c('ALPL','LMCD1','PLIN1','ADIPOQ','ENG','LEPR','THY1')
for (i in 1 : length(GOI)){
  intgroup = "Group"
  gene <- RNA_tbl[RNA_tbl$Symbol == GOI[i],'RefSeqID']
  counts <-plotCounts(dds_mrg,gene,intgroup, returnData = TRUE)
  counts$Group <- factor(counts$Group, levels=c("3dOb_4dAd_4dOb", "3dOb_4dAd", "3dOb", "Msc", "3dAd", "3dAd_4dOb","3dAd_4dOb_4dAd"))
  counts <- na.omit(counts)
  p <- ggbarplot(counts, x="Group", y="count", main = GOI[i], add = c("mean_se", "jitter"))
  print(p)
  #Extract p-values
  print(paste(GOI[i],":",RNA_tbl[RNA_tbl$Symbol==GOI[i],c("padj_3dOb_4dAd_4dOb_vs_3dOb_4dAd","padj_3dOb_4dAd_vs_3dOb","padj_3dOb_vs_Msc","padj_3dAd_vs_Msc","padj_3dAd_4dOb_vs_3dAd","padj_3dAd_4dOb_4dAd_vs_3dAd_4dOb")]))
}

#Extract p-values


