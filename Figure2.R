### Figure 2
# load the following objects from https://osf.io/c84vs/
dds_mrg <- readRDS("dds_mrg.Rds")
RNA_tbl <- readRDS("RNA_tbl3.Rds")
TiCoNE <- read.delim("TiCoNE_RefSeqID_cluster.txt")
  
### Figure 2B
# Barplots for gene of interest (GOI)

#Find Ref-seq ID's of GOI 
ALPL <- RNA_tbl[grep('ALPL\\|',RNA_tbl$Annotation.Divergence.x),'RefSeqID.x']
LMCD1 <- RNA_tbl[grep('LMCD1\\|',RNA_tbl$Annotation.Divergence.x),'RefSeqID.x']
PLIN1 <- RNA_tbl[grep('PLIN1\\|',RNA_tbl$Annotation.Divergence.x),'RefSeqID.x']
ADIPOQ <- RNA_tbl[grep('ADIPOQ\\|',RNA_tbl$Annotation.Divergence.x),'RefSeqID.x']

#Determine gene counts of GOI in each sample and do the barplot
gene <- ALPL #(Otherwise LMCD1, PLIN1 or ADIPOQ)
gene <- RNA_tbl[RNA_tbl$RefSeqID.x == gene,'RefSeqID.x']
intgroup = "Group"
counts <-plotCounts(dds_mrg,gene,intgroup, returnData = TRUE)
counts$Group <- factor(counts$Group, levels=c("7dOb_4dAd", "7dOb", "3dOb", "1dOb", "4hOb", "Msc","4hAd", "1dAd", "3dAd", "7dAd", "7dAd_4dOb"))
counts <- na.omit(counts)
ggbarplot(counts, x="Group", y="count", main = "GOI", add = c("mean_se", "jitter"))

# Final aesthetics were done in Illustrator

### Figure 2C
# PCA plot of normal diff up to day 7 and changing media ones
library(DESeq2)
library(genefilter)
library(ggplot2)

# Get rlog values of samples of interest and subject to PCA
rld_RNA <- data.frame(assay(vst(dds_mrg, blind=F)))
rld_RNA <- rld_RNA[,c(grep("_1",colnames(rld_RNA)),grep("_2",colnames(rld_RNA)))]
colnames(rld_RNA)
rld_RNA <- rld_RNA[,c(2:10,16,17,19:27,33,34)]

# PCA plot
object <- rld_RNA
intgroup <- c('Group','Replicate')
intgroup.df <- data.frame(colData(dds_mrg))
intgroup.df <- intgroup.df[colnames(object),]
ntop=300

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

### Figure 2D
# Scaled line plots for cluster specific gene expression from Rauch et al., https://doi.org/10.1038/s41588-019-0359-1 
library(DESeq2)

# Get rlog values of samples of interest
rld_RNA <- data.frame(assay(vst(dds_mrg, blind=F)))
rld_RNA <- rld_RNA[,c(grep("_1",colnames(rld_RNA)),grep("_2",colnames(rld_RNA)))]
colnames(rld_RNA)
# order samples
rld_RNA <- rld_RNA[,c(16,2:10,17,33,19:27,33,34)]
# get average
rld_RNA <- (rld_RNA[,1:11] + rld_RNA[,12:22])/2
# scale data
rld_RNA <- data.frame(t(scale(t(rld_RNA))))

#Function for transparency
addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

# plot individual genes and average for 3 selected clusters
par(mfrow=c(1,3), pty="s")
for (i in c(5,9,1)){
  genes <- TiCoNE[TiCoNE$tz_cluster_new == i,"RefSeqID"]
  genes <- genes[genes %in% rownames(rld_RNA)]
  tmp <- rld_RNA[genes,]
  Mycol <- addTrans(rainbow(14)[i],255*40 / dim(tmp)[1])
  MyGrey <- addTrans('grey',255*40 / dim(tmp)[1])
  plot(0,0,pch="", xlim=c(1,11), ylim=c(-4,4), xlab="", ylab="", yaxt="n", xaxt="n")
  mtext("Scaled RNA", side=2, line=2.5, cex=0.5)
  axis(2, at=c(-2,0,2), lab=c(-2,0,2), cex.axis=1)
  axis(1, at=c(1:11), lab=FALSE)
  for (m in 1:nrow(tmp)){
    lines(2:10,tmp[m,2:10],col=Mycol)
    lines(1:2,tmp[m,1:2], col=MyGrey)
    lines(10:11,tmp[m,10:11], col=MyGrey)
  }
  lines(1:11, apply(tmp,2, median), lwd=3)
}

# Final aesthetics were done in Illustrator