### Figure 7
# load the following objects from https://osf.io/c84vs/
dds_mrg <- readRDS("dds_mrg.Rds")
RNA_tbl <- readRDS("RNA_tbl.Rds")
# final.Confident_MED1_noadj.txt is from published enhancer dynamics of hBM-MSC-TERT4 cells (https://doi.org/10.1038/s41588-019-0359-1)
Enhancer <- read.delim("final.Confident_MED1_noadj.txt",h=T)


### Figure 7B
# Define genes regulated during adipocyte differentiation, and whehter they are less, similar, or more induced if cells were prestimulated with Ob media
# Gene groups based on 3dAd and 3dOb + 4dAd (more, less, similar induction or repression) and show them in a heatmap
library(fields)
library(pheatmap)
library(DESeq2)
library(RColorBrewer)

# Get rlog values
tmp <- data.frame(assay(rlogTransformation(dds_mrg, blind = FALSE)))

# Focus on induced Ad genes and split whether they are regulated during Ob differentiation (Common) or not (Spec) 
genes <- RNA_tbl[grepl("_Ind_",RNA_tbl$clust_3dAd),]
genes$clust_3dAd <- factor(genes$clust_3dAd, levels=c("Spec_Ind_more","Spec_Ind_similar","Spec_Ind_less","Common_Ind_more","Common_Ind_similar","Common_Ind_less"))
genes <- genes[order(genes$clust_3dAd),]

# Split induced groups into two factors, type (Spec/Common) and trend (more, similar, less)
genes$Type  <- factor(ifelse(grepl("^Spec",   genes$clust_3dAd), "Spec",  "Common"),
                      levels = c("Spec","Common"))
genes$Trend <- factor(ifelse(grepl("more",    genes$clust_3dAd), "more",
                             ifelse(grepl("less",   genes$clust_3dAd), "less", "similar")),
                      levels = c("more","similar","less"))

# Make heatmap based on genes and samples for Msc, 3dOb, 3dAd, and 3dOb + 4dAd
y2 <- t(scale(t(tmp[genes$RefSeqID,c("X3dOb_1","X3dOb_2","X3dOb_3","Msc_1","Msc_2","Msc_3","X3dAd_1","X3dAd_2","X3dAd_3","X3dOb_4dAd_1","X3dOb_4dAd_2")])))
y2[y2 > 2] <- 2
y2[y2 < -2] <- -2

# Make annotation for the clusters as a row side color
ann_row <- data.frame(
  Type  = genes$Type,
  Trend = genes$Trend
)
rownames(ann_row) <- genes$RefSeqID

#Colors for type and trend
col_trend <- c(more = "red", similar = "yellow", less = "blue")
col_type  <- c(Spec = "white", Common = "grey")

ann_colors <- list(
  Type  = col_type,
  Trend = col_trend
)

# Heatmap
mat_col <- designer.colors(n=100,col=brewer.pal(9, "Spectral"))
pheatmap(
  mat               = as.matrix(y2),
  color             = rev(mat_col),
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  treeheight_row    = 0,
  treeheight_col    = 0,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  border_color      = NA,
  annotation_row    = ann_row,
  annotation_colors = ann_colors,
)
nrow(genes[genes$Type =="Spec",])
nrow(genes[genes$Type =="Common",])

### Figure 7C
# Barplots for ABCA9 (Ad-specifically induced more) and PNMA2 (Ad-specifically induced less)
library(ggpubr)
library(DESeq2)

# Define GOI
GOI <- c('ABCA9','PNMA2')
for (i in 1 : length(GOI)){
  intgroup = "Group"
  gene <- RNA_tbl[RNA_tbl$Symbol == GOI[i],'RefSeqID']
  counts <-plotCounts(dds_mrg,gene,intgroup, returnData = TRUE)
  counts$Group <- factor(counts$Group, levels=c("3dOb_4dAd","3dOb", "Msc", "3dAd"))
  counts <- na.omit(counts)
  p <- ggbarplot(counts, x="Group", y="count", main = GOI[i], add = c("mean_se", "jitter"))
  print(p)
  #Extract p-values
  print(paste(GOI[i],":",RNA_tbl[RNA_tbl$Symbol==GOI[i],c("padj_3dOb_vs_Msc","padj_3dAd_vs_Msc","padj_3dOb_4dAd_vs_Msc","padj_3dOb_4dAd_vs_3dOb")]))
}


### Figure 7E
# Enhancer dynamics near gene groups (fold change in MED1 ChIP-seq or DNase-seq comparing 3d Ob versus Msc)

# Define induced groups, more up, similar up, less up
Spec_Ind_more <- RNA_tbl[RNA_tbl$clust_3dAd=="Spec_Ind_more","RefSeqID"]
Spec_Ind_similar <- RNA_tbl[RNA_tbl$clust_3dAd=="Spec_Ind_similar","RefSeqID"]
Spec_Ind_less <- RNA_tbl[RNA_tbl$clust_3dAd=="Spec_Ind_less","RefSeqID"]
Common_Ind_more <- RNA_tbl[RNA_tbl$clust_3dAd=="Common_Ind_more","RefSeqID"]
Common_Ind_similar <- RNA_tbl[RNA_tbl$clust_3dAd=="Common_Ind_similar","RefSeqID"]
Common_Ind_less <- RNA_tbl[RNA_tbl$clust_3dAd=="Common_Ind_less","RefSeqID"]

# Plot DNase-seq changes
par(mfrow=c(1,2), pty="s")
boxplot(
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_DHS_Ob3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"log2FC_DHS_Ob3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_DHS_Ob3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"log2FC_DHS_Ob3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_DHS_Ob3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"log2FC_DHS_Ob3d"],
  names=c("Smore","Cmore","Ssimilar","Csimilar","Sless","Cless"), col=c('white','grey'),outline=F, ylim=c(-2,3), ylab="DNase1 3dOb")
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_DHS_Ob3d"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"log2FC_DHS_Ob3d"])$estimate,
            "more"), line=2.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_DHS_Ob3d"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"log2FC_DHS_Ob3d"])$estimate,
            "similar"), line=1.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_DHS_Ob3d"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"log2FC_DHS_Ob3d"])$estimate,
            "less"), line=0.4)

effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_DHS_Ob3d"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_DHS_Ob3d"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_DHS_Ob3d"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_DHS_Ob3d"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_DHS_Ob3d"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_DHS_Ob3d"])$estimate

# Plot MED1 ChIP-seq changes
boxplot(
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_MED1_Ob3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"log2FC_MED1_Ob3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_MED1_Ob3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"log2FC_MED1_Ob3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_MED1_Ob3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"log2FC_MED1_Ob3d"],
  names=c("Smore","Cmore","Ssimilar","Csimilar","Sless","Cless"), col=c('white','grey'),outline=F, ylim=c(-2,3), ylab="Med1 3dOb")
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_MED1_Ob3d"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"log2FC_MED1_Ob3d"])$estimate,
            "more"), line=2.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_MED1_Ob3d"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"log2FC_MED1_Ob3d"])$estimate,
            "similar"), line=1.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_MED1_Ob3d"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"log2FC_MED1_Ob3d"])$estimate,
            "less"), line=0.4)

effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_MED1_Ob3d"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_MED1_Ob3d"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_MED1_Ob3d"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_MED1_Ob3d"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_MED1_Ob3d"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_MED1_Ob3d"])$estimate


### Figure 7F
# Enhancer activity near gene groups (MED1 ChIP-seq or DNase-seq in Msc)

# Define induced groups, more up, similar up, less up
Spec_Ind_more <- RNA_tbl[RNA_tbl$clust_3dAd=="Spec_Ind_more","RefSeqID"]
Spec_Ind_similar <- RNA_tbl[RNA_tbl$clust_3dAd=="Spec_Ind_similar","RefSeqID"]
Spec_Ind_less <- RNA_tbl[RNA_tbl$clust_3dAd=="Spec_Ind_less","RefSeqID"]
Common_Ind_more <- RNA_tbl[RNA_tbl$clust_3dAd=="Common_Ind_more","RefSeqID"]
Common_Ind_similar <- RNA_tbl[RNA_tbl$clust_3dAd=="Common_Ind_similar","RefSeqID"]
Common_Ind_less <- RNA_tbl[RNA_tbl$clust_3dAd=="Common_Ind_less","RefSeqID"]

# Plot DNase-seq tag counts 
par(mfrow=c(1,2), pty="s")
boxplot(
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_DHS_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"TagCount_DHS_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_DHS_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"TagCount_DHS_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_DHS_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"TagCount_DHS_Msc"],
  names=c("Smore","Cmore","Ssimilar","Csimilar","Sless","Cless"), col=c('white','grey'),outline=F,  ylab="DNase1 Msc")
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_DHS_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"TagCount_DHS_Msc"])$estimate,
            "more"), line=2.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_DHS_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"TagCount_DHS_Msc"])$estimate,
            "similar"), line=1.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_DHS_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"TagCount_DHS_Msc"])$estimate,
            "less"), line=0.4)

effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_DHS_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_DHS_Msc"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_DHS_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_DHS_Msc"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_DHS_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_DHS_Msc"])$estimate

# Plot MED1 ChIP-seq tag counts
boxplot(
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_MED1_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"TagCount_MED1_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_MED1_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"TagCount_MED1_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_MED1_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"TagCount_MED1_Msc"],
  names=c("Smore","Cmore","Ssimilar","Csimilar","Sless","Cless"), col=c('white','grey'),outline=F,  ylab="MED1 Msc")
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_MED1_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"TagCount_MED1_Msc"])$estimate,
            "more"), line=2.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_MED1_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"TagCount_MED1_Msc"])$estimate,
            "similar"), line=1.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_MED1_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"TagCount_MED1_Msc"])$estimate,
            "less"), line=0.4)

effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_MED1_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_MED1_Msc"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_MED1_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_MED1_Msc"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_MED1_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_MED1_Msc"])$estimate


### Figure 7G
# Define genes regulated during osteoblast differentiation, and whehter they are less, similar, or more induced if cells were prestimulated with Ad media
# Gene groups based on 3dOb and 3dAd + 4dOb (more, less, similar induction or repression) and show them in a heatmap
library(fields)
library(pheatmap)
library(DEseq2)
library(RColorBrewer)

# Get rlog values
tmp <- data.frame(assay(rlogTransformation(dds_mrg, blind = FALSE)))

# Focus on induced Ob genes and split whether they are regulated during Ad differentiation (Common) or not (Spec) 
genes <- RNA_tbl[grepl("_Ind_",RNA_tbl$clust_3dOb),]
genes$clust_3dOb <- factor(genes$clust_3dOb, levels=c("Spec_Ind_more","Spec_Ind_similar","Spec_Ind_less","Common_Ind_more","Common_Ind_similar","Common_Ind_less"))
genes <- genes[order(genes$clust_3dOb),]

# Split induced groups into two factors, type (Spec/Common) and trend (more, similar, less)
genes$Type  <- factor(ifelse(grepl("^Spec",   genes$clust_3dOb), "Spec",  "Common"),
                      levels = c("Spec","Common"))
genes$Trend <- factor(ifelse(grepl("more",    genes$clust_3dOb), "more",
                             ifelse(grepl("less",   genes$clust_3dOb), "less", "similar")),
                      levels = c("more","similar","less"))

# Make heatmap based on genes and samples for Msc, 3dOb, 3dAd, and 3dOb + 4dAd
y2 <- t(scale(t(tmp[genes$RefSeqID,c("X3dAd_4dOb_1","X3dAd_4dOb_2","X3dOb_1","X3dOb_2","X3dOb_3","Msc_1","Msc_2","Msc_3","X3dAd_1","X3dAd_2","X3dAd_3")])))
y2[y2 > 2] <- 2
y2[y2 < -2] <- -2

# Make annotation for the clusters as a row side color
ann_row <- data.frame(
  Type  = genes$Type,
  Trend = genes$Trend
)
rownames(ann_row) <- genes$RefSeqID

#Colors for type and trend
col_trend <- c(more = "red", similar = "yellow", less = "blue")
col_type  <- c(Spec = "white", Common = "grey")

ann_colors <- list(
  Type  = col_type,
  Trend = col_trend
)

# Heatmap
mat_col <- designer.colors(n=100,col=brewer.pal(9, "Spectral"))
pheatmap(
  mat               = as.matrix(y2),
  color             = rev(mat_col),
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  treeheight_row    = 0,
  treeheight_col    = 0,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  border_color      = NA,
  annotation_row    = ann_row,
  annotation_colors = ann_colors,
)
nrow(genes[genes$Type =="Spec",])
nrow(genes[genes$Type =="Common",])


### Figure 7H
# Barplots for ABCA9 (Ad-specifically induced more) and PNMA2 (Ad-specifically induced less)
library(ggpubr)
library(DESeq2)

# Define GOI
GOI <- c('ALPL','RANBP9')
for (i in 1 : length(GOI)){
  intgroup = "Group"
  gene <- RNA_tbl[RNA_tbl$Symbol == GOI[i],'RefSeqID']
  counts <-plotCounts(dds_mrg,gene,intgroup, returnData = TRUE)
  counts$Group <- factor(counts$Group, levels=c("3dOb", "Msc", "3dAd","3dAd_4dOb"))
  counts <- na.omit(counts)
  p <- ggbarplot(counts, x="Group", y="count", main = GOI[i], add = c("mean_se", "jitter"))
  print(p)
  #Extract p-values
  print(paste(GOI[i],":",RNA_tbl[RNA_tbl$Symbol==GOI[i],c("padj_3dOb_vs_Msc","padj_3dAd_vs_Msc","padj_3dAd_4dOb_vs_3dAd","padj_3dAd_4dOb_vs_Msc")]))
}


### Figure 7J
# Enhancer dynamics near gene groups (fold change in MED1 ChIP-seq or DNase-seq comparing 3d Ad versus Msc)

# Define induced groups, more up, similar up, less up
Spec_Ind_more <- RNA_tbl[RNA_tbl$clust_3dOb=="Spec_Ind_more","RefSeqID"]
Spec_Ind_similar <- RNA_tbl[RNA_tbl$clust_3dOb=="Spec_Ind_similar","RefSeqID"]
Spec_Ind_less <- RNA_tbl[RNA_tbl$clust_3dOb=="Spec_Ind_less","RefSeqID"]
Common_Ind_more <- RNA_tbl[RNA_tbl$clust_3dOb=="Common_Ind_more","RefSeqID"]
Common_Ind_similar <- RNA_tbl[RNA_tbl$clust_3dOb=="Common_Ind_similar","RefSeqID"]
Common_Ind_less <- RNA_tbl[RNA_tbl$clust_3dOb=="Common_Ind_less","RefSeqID"]

# Plot DNase-seq changes
par(mfrow=c(1,2), pty="s")
boxplot(
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_DHS_Ad3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"log2FC_DHS_Ad3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_DHS_Ad3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"log2FC_DHS_Ad3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_DHS_Ad3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"log2FC_DHS_Ad3d"],
  names=c("Smore","Cmore","Ssimilar","Csimilar","Sless","Cless"), col=c('white','grey'),outline=F, ylim=c(-3,5), ylab="DNase1 3dAd")
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_DHS_Ad3d"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"log2FC_DHS_Ad3d"])$estimate,
            "more"), line=2.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_DHS_Ad3d"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"log2FC_DHS_Ad3d"])$estimate,
            "similar"), line=1.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_DHS_Ad3d"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"log2FC_DHS_Ad3d"])$estimate,
            "less"), line=0.4)

effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_DHS_Ad3d"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_DHS_Ad3d"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_DHS_Ad3d"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_DHS_Ad3d"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_DHS_Ad3d"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_DHS_Ad3d"])$estimate

# Plot MED1 ChIP-seq changes
boxplot(
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_MED1_Ad3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"log2FC_MED1_Ad3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_MED1_Ad3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"log2FC_MED1_Ad3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_MED1_Ad3d"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"log2FC_MED1_Ad3d"],
  names=c("Smore","Cmore","Ssimilar","Csimilar","Sless","Cless"), col=c('white','grey'),outline=F, ylim=c(-3,5), ylab="Med1 3dAd")
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_MED1_Ad3d"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"log2FC_MED1_Ad3d"])$estimate,
            "more"), line=2.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_MED1_Ad3d"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"log2FC_MED1_Ad3d"])$estimate,
            "similar"), line=1.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_MED1_Ad3d"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"log2FC_MED1_Ad3d"])$estimate,
            "less"), line=0.4)

effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_MED1_Ad3d"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_MED1_Ad3d"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"log2FC_MED1_Ad3d"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_MED1_Ad3d"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"log2FC_MED1_Ad3d"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"log2FC_MED1_Ad3d"])$estimate


### Figure 7K
# Enhancer activity near gene groups (MED1 ChIP-seq or DNase-seq in Msc)

# Define induced groups, more up, similar up, less up
Spec_Ind_more <- RNA_tbl[RNA_tbl$clust_3dOb=="Spec_Ind_more","RefSeqID"]
Spec_Ind_similar <- RNA_tbl[RNA_tbl$clust_3dOb=="Spec_Ind_similar","RefSeqID"]
Spec_Ind_less <- RNA_tbl[RNA_tbl$clust_3dOb=="Spec_Ind_less","RefSeqID"]
Common_Ind_more <- RNA_tbl[RNA_tbl$clust_3dOb=="Common_Ind_more","RefSeqID"]
Common_Ind_similar <- RNA_tbl[RNA_tbl$clust_3dOb=="Common_Ind_similar","RefSeqID"]
Common_Ind_less <- RNA_tbl[RNA_tbl$clust_3dOb=="Common_Ind_less","RefSeqID"]

# Plot DNase-seq tag counts 
par(mfrow=c(1,2),pty="s")
boxplot(
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_DHS_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"TagCount_DHS_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_DHS_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"TagCount_DHS_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_DHS_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"TagCount_DHS_Msc"],
  names=c("Smore","Cmore","Ssimilar","Csimilar","Sless","Cless"), col=c('white','grey'),outline=F,  ylab="DNase1 Msc")
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_DHS_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"TagCount_DHS_Msc"])$estimate,
            "more"), line=2.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_DHS_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"TagCount_DHS_Msc"])$estimate,
            "similar"), line=1.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_DHS_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"TagCount_DHS_Msc"])$estimate,
            "less"), line=0.4)

effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_DHS_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_DHS_Msc"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_DHS_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_DHS_Msc"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_DHS_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_DHS_Msc"])$estimate

# Plot MED1 ChIP-seq tag counts 
boxplot(
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_MED1_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"TagCount_MED1_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_MED1_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"TagCount_MED1_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_MED1_Msc"],
  Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"TagCount_MED1_Msc"],
  names=c("Smore","Cmore","Ssimilar","Csimilar","Sless","Cless"), col=c('white','grey'),outline=F,  ylab="MED1 Msc")
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_MED1_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_more,"TagCount_MED1_Msc"])$estimate,
            "more"), line=2.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_MED1_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_similar,"TagCount_MED1_Msc"])$estimate,
            "similar"), line=1.4)
title(paste("Eff size:",
            effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_MED1_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Common_Ind_less,"TagCount_MED1_Msc"])$estimate,
            "less"), line=0.4)

effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_MED1_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_MED1_Msc"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_more,"TagCount_MED1_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_MED1_Msc"])$estimate
effsize::cohen.d(Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_similar,"TagCount_MED1_Msc"],Enhancer[Enhancer$Nearest.Refseq %in% Spec_Ind_less,"TagCount_MED1_Msc"])$estimate
