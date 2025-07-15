### Data preparation
#load libraries
library(hexbin)
library(SingleCellExperiment) 
library(simpleSingleCell) 
library(BiocSingular) 
library(umap) 
library(Seurat)
library(Matrix) 
library(valiDrops)
library(scran) 
library(scater) 
library(DropletUtils) 
library(batchelor)
library(harmony)
library(ComplexHeatmap)
library(circlize)
library(MAST)
library(limma)
library(RANN)
library(biomaRt)
library(kBET) 
library(lisi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
library(rWikiPathways)
library(GEOquery)
library(edgeR)
library(glmnet) 
library(velocyto.R) 
library(phateR)
library(ElPiGraph.R)
library(slingshot)
library(monocle3)
library(tradeSeq)
library(seriation)
library(scuttle)
library(scDblFinder)
library(SeuratWrappers)
library(devtools)
library(goseq)
library(DESeq2) 
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(DEGreport)
library(DOSE)
library(tidyverse)
library(fields)
library(ggpubr)
library(adegenet) 
library(destiny) 
library(Hmisc)
library(naniar)
library(scales)
library(monocle)
library(VennDiagram)
library(ggVennDiagram)
library(data.table)
library(ggplot2)
library(readxl)
library(phylobase)
library(ballgown)
library(SeuratData)
library(SeuratDisk)

  ### Deseq2 for bulk RNA-seq
#Importing the data (txt-file):
Raw_Diff <- read.delim("Counts_normalDiff_TERT.txt")

Raw_MeCh <- read.delim("Counts_Media_Change.txt")

#Create an empty column named Symbol in RNA_tbl
Raw_Diff$Symbol <- NA
Raw_MeCh$Symbol <- NA

#Split the existing annotation table into characters and take the first ([1]) of the characters to fill in the column Symbol
#For Raw_Diff
for (i in 1:nrow(Raw_Diff)) 
  
{ 
  Raw_Diff[i, "Symbol"] <- unlist(strsplit(as.character(Raw_Diff[i, "Annotation.Divergence"]), "\\|"))[1]
}
#For Raw_MeCh
for (i in 1:nrow(Raw_MeCh)) 
  
{ 
  Raw_MeCh[i, "Symbol"] <- unlist(strsplit(as.character(Raw_MeCh[i, "Annotation.Divergence"]), "\\|"))[1]
}

#Deleting genes with values equal to or less than 2
Raw_Diff <- Raw_Diff[apply(Raw_Diff[,9:41],1,max)> 2,]
Raw_MeCh <- Raw_MeCh[apply(Raw_MeCh[,9:20],1,max)> 2,]

#Merging files using Symbol instead of RefSeqID
Raw_mrg <- merge(Raw_Diff, Raw_MeCh[,9:20], by = "Symbol")

#Running the DESeq
row.names(Raw_mrg) <- Raw_mrg$RefSeqID.x
df_cnt <- as.data.frame(Raw_mrg[,c(10:42,51:62)])
df_cl <- read.delim("colData_merged.txt")
dds_mrg <- DESeqDataSetFromMatrix(countData = df_cnt, colData = df_cl, design = ~ Group)
dds_mrg <- DESeq(dds_mrg)
#saveRDS(dds_mrg, "dds_mrg.Rds")

#Create a data table showing Pval, Padj & log2FoldChange of chosen comparisons
RNA_tbl <- Raw_mrg[, c(1:7, 9)]
RNA_tbl <- cbind(RNA_tbl, 
                 as.data.frame(results(dds_mrg, contrast = c("Group", "14dOb", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dOb", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "1dOb", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "4hOb", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "4hAd", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "1dAd", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dAd", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "14dAd", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb_4dAd", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd_4dOb", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb_4dAd_4dOb", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd_4dOb_4dAd", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dOb_4dAd", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dAd_4dOb", "Msc")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dOb_4dAd", "7dOb")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dAd_4dOb", "7dAd")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dOb_4dAd", "3dAd")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd_4dOb", "3dAd")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb_4dAd", "3dOb")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd_4dOb_4dAd", "3dAd_4dOb")))[5],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb_4dAd_4dOb", "3dOb_4dAd")))[5],
                 
                 as.data.frame(results(dds_mrg, contrast = c("Group", "14dOb", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dOb", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "1dOb", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "4hOb", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "4hAd", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "1dAd", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dAd", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "14dAd", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb_4dAd", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd_4dOb", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb_4dAd_4dOb", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd_4dOb_4dAd", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dOb_4dAd", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dAd_4dOb", "Msc")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dOb_4dAd", "7dOb")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dAd_4dOb", "7dAd")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd_4dOb", "3dAd")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb_4dAd", "3dOb")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd_4dOb_4dAd", "3dAd_4dOb")))[6],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb_4dAd_4dOb", "3dOb_4dAd")))[6],
                 
                 as.data.frame(results(dds_mrg, contrast = c("Group", "14dOb", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dOb", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "1dOb", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "4hOb", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "4hAd", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "1dAd", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dAd", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "14dAd", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb_4dAd", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd_4dOb", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb_4dAd_4dOb", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd_4dOb_4dAd", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dOb_4dAd", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dAd_4dOb", "Msc")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dOb_4dAd", "7dOb")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "7dAd_4dOb", "7dAd")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd_4dOb", "3dAd")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb_4dAd", "3dOb")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dAd_4dOb_4dAd", "3dAd_4dOb")))[2],
                 as.data.frame(results(dds_mrg, contrast = c("Group", "3dOb_4dAd_4dOb", "3dOb_4dAd")))[2])

#Naming the columns of comparison
names(RNA_tbl)[9:65] <- c("pval_14dOb_vs_Msc","pval_7dOb_vs_Msc","pval_3dOb_vs_Msc","pval_1dOb_vs_Msc","pval_4hOb_vs_Msc",
                          "pval_4hAd_vs_Msc", "pval_1dAd_vs_Msc","pval_3dAd_vs_Msc","pval_7dAd_vs_Msc","pval_14dAd_vs_Msc",
                          "pval_3dOb_4dAd_vs_Msc","pval_3dAd_4dOb_vs_Msc","pval_3dOb_4dAd_4dOb_vs_Msc","pval_3dAd_4dOb_4dAd_vs_Msc",
                          "pval_7dOb_4dAd_vs_Msc","pval_7dAd_4dOb_vs_Msc", "pval_7dOb_4dAd_vs_7dOb","pval_7dAd_4dOb_vs_7dAd","pval_3dAd_4dOb_vs_3dAd","pval_3dOb_4dAd_vs_3dOb","pval_3dAd_4dOb_4dAd_vs_3dAd_4dOb","pval_3dOb_4dAd_4dOb_vs_3dOb_4dAd",
                          "padj_14dOb_vs_Msc","padj_7dOb_vs_Msc","padj_3dOb_vs_Msc","padj_1dOb_vs_Msc","padj_4hOb_vs_Msc",
                          "padj_4hAd_vs_Msc", "padj_1dAd_vs_Msc","padj_3dAd_vs_Msc","padj_7dAd_vs_Msc","padj_14dAd_vs_Msc",
                          "padj_3dOb_4dAd_vs_Msc","padj_3dAd_4dOb_vs_Msc","padj_3dOb_4dAd_4dOb_vs_Msc","padj_3dAd_4dOb_4dAd_vs_Msc",
                          "padj_7dOb_4dAd_vs_Msc","padj_7dAd_4dOb_vs_Msc", "padj_7dOb_4dAd_vs_7dOb","padj_7dAd_4dOb_vs_7dAd","padj_3dAd_4dOb_vs_3dAd","padj_3dOb_4dAd_vs_3dOb","padj_3dAd_4dOb_4dAd_vs_3dAd_4dOb","padj_3dOb_4dAd_4dOb_vs_3dOb_4dAd",
                          "logFC_14dOb_vs_Msc","logFC_7dOb_vs_Msc","logFC_3dOb_vs_Msc","logFC_1dOb_vs_Msc","logFC_4hOb_vs_Msc",
                          "logFC_4hAd_vs_Msc", "logFC_1dAd_vs_Msc","logFC_3dAd_vs_Msc","logFC_7dAd_vs_Msc","logFC_14dAd_vs_Msc",
                          "logFC_3dOb_4dAd_vs_Msc","logFC_3dAd_4dOb_vs_Msc","logFC_3dOb_4dAd_4dOb_vs_Msc","logFC_3dAd_4dOb_4dAd_vs_Msc",
                          "logFC_7dOb_4dAd_vs_Msc","logFC_7dAd_4dOb_vs_Msc", "logFC_7dOb_4dAd_vs_7dOb","logFC_7dAd_4dOb_vs_7dAd","logFC_3dAd_4dOb_vs_3dAd","logFC_3dOb_4dAd_vs_3dOb","logFC_3dAd_4dOb_4dAd_vs_3dAd_4dOb","logFC_3dOb_4dAd_4dOb_vs_3dOb_4dAd")

#Changing NA values into 1
RNA_tbl[!complete.cases(RNA_tbl$pval_14dOb_vs_Msc), "pval_14dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_7dOb_vs_Msc), "pval_7dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_3dOb_vs_Msc), "pval_3dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_1dOb_vs_Msc), "pval_1dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_4hOb_vs_Msc), "pval_4hOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_4hAd_vs_Msc), "pval_4hAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_1dAd_vs_Msc), "pval_1dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_3dAd_vs_Msc), "pval_3dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_7dAd_vs_Msc), "pval_7dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_14dAd_vs_Msc), "pval_14dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_3dOb_4dAd_vs_Msc), "pval_3dOb_4dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_3dAd_4dOb_vs_Msc), "pval_3dAd_4dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_3dOb_4dAd_4dOb_vs_Msc), "pval_3dOb_4dAd_4dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_3dAd_4dOb_4dAd_vs_Msc), "pval_3dAd_4dOb_4dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_7dOb_4dAd_vs_Msc), "pval_7dOb_4dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_7dAd_4dOb_vs_Msc), "pval_7dAd_4dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_7dOb_4dAd_vs_7dOb), "pval_7dOb_4dAd_vs_7dOb"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_7dAd_4dOb_vs_7dAd), "pval_7dAd_4dOb_vs_7dAd"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_3dAd_4dOb_vs_3dAd), "pval_3dAd_4dOb_vs_3dAd"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_3dOb_4dAd_vs_3dOb), "pval_3dOb_4dAd_vs_3dOb"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_3dAd_4dOb_4dAd_vs_3dAd_4dOb), "pval_3dAd_4dOb_4dAd_vs_3dAd_4dOb"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$pval_3dOb_4dAd_4dOb_vs_3dOb_4dAd), "pval_3dOb_4dAd_4dOb_vs_3dOb_4dAd"] <- 1

RNA_tbl[!complete.cases(RNA_tbl$padj_14dOb_vs_Msc), "padj_14dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_7dOb_vs_Msc), "padj_7dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_3dOb_vs_Msc), "padj_3dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_1dOb_vs_Msc), "padj_1dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_4hOb_vs_Msc), "padj_4hOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_4hAd_vs_Msc), "padj_4hAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_1dAd_vs_Msc), "padj_1dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_3dAd_vs_Msc), "padj_3dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_7dAd_vs_Msc), "padj_7dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_14dAd_vs_Msc), "padj_14dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_3dOb_4dAd_vs_Msc), "padj_3dOb_4dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_3dAd_4dOb_vs_Msc), "padj_3dAd_4dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_3dOb_4dAd_4dOb_vs_Msc), "padj_3dOb_4dAd_4dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_3dAd_4dOb_4dAd_vs_Msc), "padj_3dAd_4dOb_4dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_7dOb_4dAd_vs_Msc), "padj_7dOb_4dAd_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_7dAd_4dOb_vs_Msc), "padj_7dAd_4dOb_vs_Msc"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_7dOb_4dAd_vs_7dOb), "padj_7dOb_4dAd_vs_7dOb"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_7dAd_4dOb_vs_7dAd), "padj_7dAd_4dOb_vs_7dAd"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_3dAd_4dOb_vs_3dAd), "padj_3dAd_4dOb_vs_3dAd"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_3dOb_4dAd_vs_3dOb), "padj_3dOb_4dAd_vs_3dOb"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_3dAd_4dOb_4dAd_vs_3dAd_4dOb), "padj_3dAd_4dOb_4dAd_vs_3dAd_4dOb"] <- 1
RNA_tbl[!complete.cases(RNA_tbl$padj_3dOb_4dAd_4dOb_vs_3dOb_4dAd), "padj_3dOb_4dAd_4dOb_vs_3dOb_4dAd"] <- 1
#saveRDS(RNA_tbl, "RNA_tbl.Rds")

  ### scRNA RNA seq seurat object
## Import count objects generated using zUMIs
TERT_MSC <- readRDS("TERT_MSC.dgecounts.rds")
TERT_7dOb <- readRDS("TERT_7dOb.dgecounts.rds")
TERT_7dAd <- readRDS("TERT_7dAd.dgecounts.rds")
TERT_7dOb_4dAd <- readRDS("TERT_7dOb_4dAd.dgecounts.rds")
TERT_7dAd_4dOb <- readRDS("TERT_7dAd_4dOb.dgecounts.rds")

# Extracting intron+exons counts 
TERT_MSC <- TERT_MSC$umicount$inex$all
TERT_7dOb <- TERT_7dOb$umicount$inex$all
TERT_7dAd <- TERT_7dAd$umicount$inex$all
TERT_7dOb_4dAd <- TERT_7dOb_4dAd$umicount$inex$all
TERT_7dAd_4dOb <- TERT_7dAd_4dOb$umicount$inex$all

## Fill in the matrices to give them all the same dimensions
Genes <- unique(c(rownames(TERT_MSC),rownames(TERT_7dOb),rownames(TERT_7dAd),rownames(TERT_7dAd_4dOb),rownames(TERT_7dOb_4dAd)))

# Insert empty lines with missing genes to get same dimensions
Tmp <- as.sparse(matrix(ncol=ncol(TERT_MSC), nrow=length(Genes[!(Genes %in% rownames(TERT_MSC))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(TERT_MSC))]
TERT_MSC <- Matrix::rbind2(TERT_MSC, Tmp)
TERT_MSC <- TERT_MSC[ order(rownames(TERT_MSC)),]
dim(TERT_MSC) # 43610 2597

Tmp <- as.sparse(matrix(ncol=ncol(TERT_7dOb), nrow=length(Genes[!(Genes %in% rownames(TERT_7dOb))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(TERT_7dOb))]
TERT_7dOb <- Matrix::rbind2(TERT_7dOb, Tmp)
TERT_7dOb <- TERT_7dOb[ order(rownames(TERT_7dOb)),]
dim(TERT_7dOb) # 43610 4434

Tmp <- as.sparse(matrix(ncol=ncol(TERT_7dAd), nrow=length(Genes[!(Genes %in% rownames(TERT_7dAd))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(TERT_7dAd))]
TERT_7dAd <- Matrix::rbind2(TERT_7dAd, Tmp)
TERT_7dAd <- TERT_7dAd[ order(rownames(TERT_7dAd)),]
dim(TERT_7dAd) # 43610 8630

Tmp <- as.sparse(matrix(ncol=ncol(TERT_7dOb_4dAd), nrow=length(Genes[!(Genes %in% rownames(TERT_7dOb_4dAd))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(TERT_7dOb_4dAd))]
TERT_7dOb_4dAd <- Matrix::rbind2(TERT_7dOb_4dAd, Tmp)
TERT_7dOb_4dAd <- TERT_7dOb_4dAd[ order(rownames(TERT_7dOb_4dAd)),]
dim(TERT_7dOb_4dAd) # 43610 11328

Tmp <- as.sparse(matrix(ncol=ncol(TERT_7dAd_4dOb), nrow=length(Genes[!(Genes %in% rownames(TERT_7dAd_4dOb))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(TERT_7dAd_4dOb))]
TERT_7dAd_4dOb <- Matrix::rbind2(TERT_7dAd_4dOb, Tmp)
TERT_7dAd_4dOb <- TERT_7dAd_4dOb[ order(rownames(TERT_7dAd_4dOb)),]
dim(TERT_7dAd_4dOb) # 43610 7251

## Paste in the experiment name into the column names to make barcodes/cell IDs unique.
colnames(TERT_MSC) <- paste("TERT_MSC",colnames(TERT_MSC), sep="")
colnames(TERT_7dOb) <- paste("TERT_7dOb",colnames(TERT_7dOb), sep="")
colnames(TERT_7dAd) <- paste("TERT_7dAd",colnames(TERT_7dAd), sep="")
colnames(TERT_7dOb_4dAd) <- paste("TERT_7dOb_4dAd",colnames(TERT_7dOb_4dAd), sep="")
colnames(TERT_7dAd_4dOb) <- paste("TERT_7dAd_4dOb",colnames(TERT_7dAd_4dOb), sep="")

## Convert Ensemble IDs to Gene Symbols
# Process gene list (downloaded from BioMart), keep only non-empty gene symbols and deduplicate.
Genes <- read.delim("mart_export_human_Ens_GRCh38.txt")
Genes <- Genes[ Genes$Gene.stable.ID %in% rownames(TERT_MSC),]
Genes <- Genes[ Genes$Gene.name !=  "",]
Genes <- Genes[ duplicated(Genes$Gene.stable.ID) == F,]
Genes <- Genes[ duplicated(Genes$Gene.name) == F,]
Genes <- Genes[ order(Genes$Gene.stable.ID),]
dim(Genes)

# Subset count matrices and replace Ensemble IDs with gene symbols
TERT_MSC <- TERT_MSC[ rownames(TERT_MSC) %in% Genes$Gene.stable.ID,]
TERT_MSC <- TERT_MSC[ order(rownames(TERT_MSC)),]
rownames(TERT_MSC) <- as.character(Genes$Gene.name)
dim(TERT_MSC) # 41500 2597

TERT_7dOb <- TERT_7dOb[ rownames(TERT_7dOb) %in% Genes$Gene.stable.ID,]
TERT_7dOb <- TERT_7dOb[ order(rownames(TERT_7dOb)),]
rownames(TERT_7dOb) <- as.character(Genes$Gene.name)
dim(TERT_7dOb) # 41500 4434

TERT_7dAd <- TERT_7dAd[ rownames(TERT_7dAd) %in% Genes$Gene.stable.ID,]
TERT_7dAd <- TERT_7dAd[ order(rownames(TERT_7dAd)),]
rownames(TERT_7dAd) <- as.character(Genes$Gene.name)
dim(TERT_7dAd) # 41500 8630

TERT_7dOb_4dAd <- TERT_7dOb_4dAd[ rownames(TERT_7dOb_4dAd) %in% Genes$Gene.stable.ID,]
TERT_7dOb_4dAd <- TERT_7dOb_4dAd[ order(rownames(TERT_7dOb_4dAd)),]
rownames(TERT_7dOb_4dAd) <- as.character(Genes$Gene.name)
dim(TERT_7dOb_4dAd) # 41500 11328

TERT_7dAd_4dOb <- TERT_7dAd_4dOb[ rownames(TERT_7dAd_4dOb) %in% Genes$Gene.stable.ID,]
TERT_7dAd_4dOb <- TERT_7dAd_4dOb[ order(rownames(TERT_7dAd_4dOb)),]
rownames(TERT_7dAd_4dOb) <- as.character(Genes$Gene.name)
dim(TERT_7dAd_4dOb) # 41500 7251

## Create SingleCellExperiment objects to use DropUtils, scran, scater, etc.
TERT_MSC_sce <- SingleCellExperiment(list(counts=TERT_MSC))
TERT_7dOb_sce <- SingleCellExperiment(list(counts=TERT_7dOb))
TERT_7dAd_sce <- SingleCellExperiment(list(counts=TERT_7dAd))
TERT_7dOb_4dAd_sce <- SingleCellExperiment(list(counts=TERT_7dOb_4dAd))
TERT_7dAd_4dOb_sce <- SingleCellExperiment(list(counts=TERT_7dAd_4dOb))

##  Observational plots
br.out <- barcodeRanks(TERT_MSC)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

br.out <- barcodeRanks(TERT_7dOb)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

br.out <- barcodeRanks(TERT_7dAd)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

br.out <- barcodeRanks(TERT_7dOb_4dAd)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

br.out <- barcodeRanks(TERT_7dAd_4dOb)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
       legend=c("knee", "inflection"))

par(mfrow=c(2,3), pty="s")
tmp1 <- apply(TERT_MSC,2,sum)
tmp2 <- colSums(TERT_MSC > 0)
plot(tmp1,tmp2,main="TERT MSC")

tmp1 <- apply(TERT_7dOb,2,sum)
tmp2 <- colSums(TERT_7dOb > 0)
plot(tmp1,tmp2,main="TERT 7dOb")

tmp1 <- apply(TERT_7dAd,2,sum)
tmp2 <- colSums(TERT_7dAd > 0)
plot(tmp1,tmp2,main="TERT 7dAd")

tmp1 <- apply(TERT_7dOb_4dAd,2,sum)
tmp2 <- colSums(TERT_7dOb_4dAd > 0)
plot(tmp1,tmp2,main="TERT 7dOb 4dAd")

tmp1 <- apply(TERT_7dAd_4dOb,2,sum)
tmp2 <- colSums(TERT_7dAd_4dOb > 0)
plot(tmp1,tmp2,main="TERT 7dAd 4dOb")

## Lower threshold were set based on being in the range of knee points and adjusted to match the expectation of recovering approximately 10.000 non-empty droplets, based on the loading of the 10X chip.
par(mfrow=c(2,3), pty="s")

e.out <- emptyDrops(counts(TERT_MSC_sce))
is.cell <- e.out$FDR <= 0.01
table(is.cell)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability",main="TERT MSC")

e.out <- emptyDrops(counts(TERT_7dOb_sce))
is.cell <- e.out$FDR <= 0.01
table(is.cell)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability",main="TERT 7dOb")

e.out <- emptyDrops(counts(TERT_7dAd_sce))
is.cell <- e.out$FDR <= 0.01
table(is.cell)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability",main="TERT 7dAd")

e.out <- emptyDrops(counts(TERT_7dOb_4dAd_sce))
is.cell <- e.out$FDR <= 0.01
table(is.cell)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability",main="TERT 7dOb_4dAd")

e.out <- emptyDrops(counts(TERT_7dAd_4dOb_sce))
is.cell <- e.out$FDR <= 0.01
table(is.cell)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability",main="TERT 7dAd_4dOb")

e.out <- emptyDrops(counts(TERT_MSC_sce))
TERT_MSC_sce <- TERT_MSC_sce[,which(e.out$FDR <= 0.01)]
e.out <- emptyDrops(counts(TERT_7dOb_sce))
TERT_7dOb_sce <- TERT_7dOb_sce[,which(e.out$FDR <= 0.01)]
e.out <- emptyDrops(counts(TERT_7dAd_sce))
TERT_7dAd_sce <- TERT_7dAd_sce[,which(e.out$FDR <= 0.01)]
e.out <- emptyDrops(counts(TERT_7dOb_4dAd_sce))
TERT_7dOb_4dAd_sce <- TERT_7dOb_4dAd_sce[,which(e.out$FDR <= 0.01)]
e.out <- emptyDrops(counts(TERT_7dAd_4dOb_sce))
TERT_7dAd_4dOb_sce <- TERT_7dAd_4dOb_sce[,which(e.out$FDR <= 0.01)]

rm(tmp1,tmp2,br.out,Tmp,Genes,o,e.out,is.cell)

## Calculate QC parameters (throws a warning, that can be ignored)
TERT_MSC_sce <- addPerCellQC(TERT_MSC_sce, subsets=list(Mito=grep(pattern = "^MT-", x = rownames(TERT_MSC_sce), value = FALSE)))
TERT_7dOb_sce <- addPerCellQC(TERT_7dOb_sce, subsets=list(Mito=grep(pattern = "^MT-", x = rownames(TERT_7dOb_sce), value = FALSE)))
TERT_7dAd_sce <- addPerCellQC(TERT_7dAd_sce, subsets=list(Mito=grep(pattern = "^MT-", x = rownames(TERT_7dAd_sce), value = FALSE)))
TERT_7dOb_4dAd_sce <- addPerCellQC(TERT_7dOb_4dAd_sce, subsets=list(Mito=grep(pattern = "^MT-", x = rownames(TERT_7dOb_4dAd_sce), value = FALSE)))
TERT_7dAd_4dOb_sce <- addPerCellQC(TERT_7dAd_4dOb_sce, subsets=list(Mito=grep(pattern = "^MT", x = rownames(TERT_7dAd_4dOb_sce), value = FALSE)))

### Histograms of quality measures for each dataset
par(mfcol=c(5,4),pty="s")
hist(TERT_MSC_sce$subsets_Mito_percent,main="TERT MSC")
hist(TERT_7dOb_sce$subsets_Mito_percent,main="TERT 7dOb")
hist(TERT_7dAd_sce$subsets_Mito_percent,main="TERT 7dAd")
hist(TERT_7dOb_4dAd_sce$subsets_Mito_percent,main="TERT 7dOb_4dAd")
hist(TERT_7dAd_4dOb_sce$subsets_Mito_percent,main="TERT 7dAd_4dOb")

hist(TERT_MSC_sce$sum,main="TERT MSC")
hist(TERT_7dOb_sce$sum,main="TERT 7dOb")
hist(TERT_7dAd_sce$sum,main="TERT 7dAd")
hist(TERT_7dOb_4dAd_sce$sum,main="TERT 7dOb_4dAd")
hist(TERT_7dAd_4dOb_sce$sum,main="TERT 7dAd_4dOb")

hist(TERT_MSC_sce$detected,main="TERT MSC")
hist(TERT_7dOb_sce$detected,main="TERT 7dOb")
hist(TERT_7dAd_sce$detected,main="TERT 7dAd")
hist(TERT_7dOb_4dAd_sce$detected,main="TERT 7dOb_4dAd")
hist(TERT_7dAd_4dOb_sce$detected,main="TERT 7dAd_4dOb")

hist(TERT_MSC_sce$sum/TERT_MSC_sce$detected,main="TERT MSC")
hist(TERT_7dOb_sce$sum/TERT_7dOb_sce$detected,main="TERT 7dOb")
hist(TERT_7dAd_sce$sum/TERT_7dAd_sce$detected,main="TERT 7dAd")
hist(TERT_7dOb_4dAd_sce$sum/TERT_7dOb_4dAd_sce$detected,main="TERT 7dOb_4dAd")
hist(TERT_7dAd_4dOb_sce$sum/TERT_7dAd_4dOb_sce$detected,main="TERT 7dAd_4dOb")

#Inspect with boxplots
box_subsets_Mito_percent <- list()
box_subsets_Mito_percent[[1]] <- TERT_MSC_sce$subsets_Mito_percent
box_subsets_Mito_percent[[2]] <- TERT_7dOb_sce$subsets_Mito_percent
box_subsets_Mito_percent[[3]] <- TERT_7dAd_sce$subsets_Mito_percent
box_subsets_Mito_percent[[4]] <- TERT_7dOb_4dAd_sce$subsets_Mito_percent
box_subsets_Mito_percent[[5]] <- TERT_7dAd_4dOb_sce$subsets_Mito_percent

box_sum <- list()
box_sum[[1]] <- TERT_MSC_sce$sum
box_sum[[2]] <- TERT_7dOb_sce$sum
box_sum[[3]] <- TERT_7dAd_sce$sum
box_sum[[4]] <- TERT_7dOb_4dAd_sce$sum
box_sum[[5]] <- TERT_7dAd_4dOb_sce$sum

box_detected <- list()
box_detected[[1]] <- TERT_MSC_sce$detected
box_detected[[2]] <- TERT_7dOb_sce$detected
box_detected[[3]] <- TERT_7dAd_sce$detected
box_detected[[4]] <- TERT_7dOb_4dAd_sce$detected
box_detected[[5]] <- TERT_7dAd_4dOb_sce$detected

box_sum_detected <- list()
box_sum_detected[[1]] <- TERT_MSC_sce$sum/TERT_MSC_sce$detected
box_sum_detected[[2]] <- TERT_7dOb_sce$sum/TERT_7dOb_sce$detected
box_sum_detected[[3]] <- TERT_7dAd_sce$sum/TERT_7dAd_sce$detected
box_sum_detected[[4]] <- TERT_7dOb_4dAd_sce$sum/TERT_7dOb_4dAd_sce$detected
box_sum_detected[[5]] <- TERT_7dAd_4dOb_sce$sum/TERT_7dAd_4dOb_sce$detected

par(mfcol=c(2,2),pty="s")
boxplot(box_subsets_Mito_percent, las=2, outline = TRUE, names = c("MSC", "7dOb", "7dAd", "7dOb_4dAd", "7dAd_4dOb"), xlab='libraries', ylab='Mito_percent', main= "Mito_percent")
boxplot(box_sum, las=2, outline = TRUE, names = c("MSC", "7dOb", "7dAd", "7dOb_4dAd", "7dAd_4dOb"), xlab='libraries', ylab='sum', main= "sum")
boxplot(box_detected, las=2, outline = TRUE, names = c("MSC", "7dOb", "7dAd", "7dOb_4dAd", "7dAd_4dOb"), xlab='libraries', ylab='detected', main= "detected")
boxplot(box_sum_detected, las=2, outline = TRUE, names = c("MSC", "7dOb", "7dAd", "7dOb_4dAd", "7dAd_4dOb"), xlab='libraries', ylab='sum / detected', main= "sum / detected")

## Threshold filtering of droplets in each dataset
## REMOVE: Droplets with more than 10% mitochondrial reads, less than 1000 UMIs, less than 500 genes or extremely high ratio between counts and genes (low complexity))
TERT_MSC_sce <- TERT_MSC_sce[,!(TERT_MSC_sce$subsets_Mito_percent > 10)] 
TERT_MSC_sce <- TERT_MSC_sce[,!(TERT_MSC_sce$sum/TERT_MSC_sce$detected > 6)] 
TERT_MSC_sce <- TERT_MSC_sce[,(TERT_MSC_sce$sum >= 200 & TERT_MSC_sce$detected >= 200)] 
dim(TERT_MSC_sce) #41500 2356

TERT_7dOb_sce <- TERT_7dOb_sce[,!(TERT_7dOb_sce$subsets_Mito_percent > 10)] 
TERT_7dOb_sce <- TERT_7dOb_sce[,!(TERT_7dOb_sce$sum/TERT_7dOb_sce$detected > 6)] 
TERT_7dOb_sce <- TERT_7dOb_sce[,(TERT_7dOb_sce$sum >= 200 & TERT_7dOb_sce$detected >= 200)] 
dim(TERT_7dOb_sce) #41500 4132

TERT_7dAd_sce <- TERT_7dAd_sce[,!(TERT_7dAd_sce$subsets_Mito_percent > 10)] 
TERT_7dAd_sce <- TERT_7dAd_sce[,!(TERT_7dAd_sce$sum/TERT_7dAd_sce$detected > 6)] 
TERT_7dAd_sce <- TERT_7dAd_sce[,(TERT_7dAd_sce$sum >= 200 & TERT_7dAd_sce$detected >= 200)] 
dim(TERT_7dAd_sce) #41500 7889

TERT_7dOb_4dAd_sce <- TERT_7dOb_4dAd_sce[,!(TERT_7dOb_4dAd_sce$subsets_Mito_percent > 10)] 
TERT_7dOb_4dAd_sce <- TERT_7dOb_4dAd_sce[,!(TERT_7dOb_4dAd_sce$sum/TERT_7dOb_4dAd_sce$detected > 6)] 
TERT_7dOb_4dAd_sce <- TERT_7dOb_4dAd_sce[,(TERT_7dOb_4dAd_sce$sum >= 200 & TERT_7dOb_4dAd_sce$detected >= 200)] 
dim(TERT_7dOb_4dAd_sce) #41500 10827

TERT_7dAd_4dOb_sce <- TERT_7dAd_4dOb_sce[,!(TERT_7dAd_4dOb_sce$subsets_Mito_percent > 10)] 
TERT_7dAd_4dOb_sce <- TERT_7dAd_4dOb_sce[,!(TERT_7dAd_4dOb_sce$sum/TERT_7dAd_4dOb_sce$detected > 6)] 
TERT_7dAd_4dOb_sce <- TERT_7dAd_4dOb_sce[,(TERT_7dAd_4dOb_sce$sum >= 200 & TERT_7dAd_4dOb_sce$detected >= 200)] 
dim(TERT_7dAd_4dOb_sce) #41500 6450

## Automatic filtering droplets in each dataset using PCA across all QC metrics
# Calculate outliers
TERT_MSC_sce <- logNormCounts(TERT_MSC_sce)
TERT_7dOb_sce <- logNormCounts(TERT_7dOb_sce)
TERT_7dAd_sce <- logNormCounts(TERT_7dAd_sce)
TERT_7dOb_4dAd_sce <- logNormCounts(TERT_7dOb_4dAd_sce)
TERT_7dAd_4dOb_sce <- logNormCounts(TERT_7dAd_4dOb_sce)

TERT_MSC_sce <- runColDataPCA(TERT_MSC_sce,variables=list("sum","detected","subsets_Mito_percent"),outliers=T)
TERT_7dOb_sce <- runColDataPCA(TERT_7dOb_sce,variables=list("sum","detected","subsets_Mito_percent"),outliers=T)
TERT_7dAd_sce <- runColDataPCA(TERT_7dAd_sce,variables=list("sum","detected","subsets_Mito_percent"),outliers=T)
TERT_7dOb_4dAd_sce <- runColDataPCA(TERT_7dOb_4dAd_sce,variables=list("sum","detected","subsets_Mito_percent"),outliers=T)
TERT_7dAd_4dOb_sce <- runColDataPCA(TERT_7dAd_4dOb_sce,variables=list("sum","detected","subsets_Mito_percent"),outliers=T)

table(TERT_MSC_sce$outlier)
table(TERT_7dOb_sce$outlier)
table(TERT_7dAd_sce$outlier)
table(TERT_7dOb_4dAd_sce$outlier)
table(TERT_7dAd_4dOb_sce$outlier)

# Remove outliers with PCA
TERT_MSC_sce <- TERT_MSC_sce[ ,!TERT_MSC_sce$outlier] 
TERT_7dOb_sce <- TERT_7dOb_sce[ ,!TERT_7dOb_sce$outlier] 
TERT_7dAd_sce <- TERT_7dAd_sce[ ,!TERT_7dAd_sce$outlier] 
TERT_7dOb_4dAd_sce <- TERT_7dOb_4dAd_sce[ ,!TERT_7dOb_4dAd_sce$outlier] 
TERT_7dAd_4dOb_sce <- TERT_7dAd_4dOb_sce[ ,!TERT_7dAd_4dOb_sce$outlier] 

## Threshold filtering of genes in each dataset
## REMOVE: Genes expressed in less than 10 nuclei in all datasets
# Find lowly expressed genes and get the intersection
TERT_MSC_low <- names(which(nexprs(TERT_MSC_sce, byrow=T) <= 1))
TERT_7dOb_low <- names(which(nexprs(TERT_7dOb_sce, byrow=T) <= 1))
TERT_7dAd_low <- names(which(nexprs(TERT_7dAd_sce, byrow=T) <= 1))
TERT_7dOb_4dAd_low <- names(which(nexprs(TERT_7dOb_4dAd_sce, byrow=T) <= 1))
TERT_7dAd_4dOb_low <- names(which(nexprs(TERT_7dAd_4dOb_sce, byrow=T) <= 1))
Low <- Reduce(intersect, list(TERT_MSC_low,TERT_7dOb_low,TERT_7dAd_low, TERT_7dOb_4dAd_low, TERT_7dAd_4dOb_low))

# Remove lowly expressed genes
TERT_MSC_sce <- TERT_MSC_sce[ which(!(rownames(TERT_MSC_sce) %in% Low)),] 
TERT_7dOb_sce <- TERT_7dOb_sce[ which(!(rownames(TERT_7dOb_sce) %in% Low)),] 
TERT_7dAd_sce <- TERT_7dAd_sce[ which(!(rownames(TERT_7dAd_sce) %in% Low)),] 
TERT_7dOb_4dAd_sce <- TERT_7dOb_4dAd_sce[ which(!(rownames(TERT_7dOb_4dAd_sce) %in% Low)),] 
TERT_7dAd_4dOb_sce <- TERT_7dAd_4dOb_sce[ which(!(rownames(TERT_7dAd_4dOb_sce) %in% Low)),] 

## Filtering genes based on biotype and transcript level support (not done yet)
## REMOVE: Non-protein coding genes, keep the one with high transcript level support (tsl)
Genes <- read.delim("mart_export_human_Ens_GRCh38.txt")
Genes <- Genes[ Genes$Gene.name %in% rownames(TERT_MSC_sce),]
Genes <- Genes[ grep("tsl1", Genes$Transcript.support.level..TSL.),]
Genes <- Genes[ Genes$Gene.name !=  "",]
Genes <- Genes[ Genes$Transcript.type == "protein_coding",]
Genes <- Genes[ Genes$Gene.type == "protein_coding",]
Genes <- Genes[ !is.na(Genes$Gene.name),]
Genes <- Genes[ duplicated(Genes$Gene.stable.ID) == F,] 

# Filter the sce objects and ambient genes
TERT_MSC_sce <- TERT_MSC_sce[ which(rownames(TERT_MSC_sce) %in% Genes$Gene.name),] # 15743 genes left
TERT_7dOb_sce <- TERT_7dOb_sce[ which(rownames(TERT_7dOb_sce) %in% Genes$Gene.name),] # 15743 genes left
TERT_7dAd_sce <- TERT_7dAd_sce[ which(rownames(TERT_7dAd_sce) %in% Genes$Gene.name),] # 15743 genes left
TERT_7dOb_4dAd_sce <- TERT_7dOb_4dAd_sce[ which(rownames(TERT_7dOb_4dAd_sce) %in% Genes$Gene.name),] # 15743 genes left
TERT_7dAd_4dOb_sce <- TERT_7dAd_4dOb_sce[ which(rownames(TERT_7dAd_4dOb_sce) %in% Genes$Gene.name),] # 15743 genes left

## Normalize the count matrices #
# Cluster each data. NOTE: THIS STEP IS NON-DETERMINISTIC - RESULTS VARY FROM RUN TO RUN
TERT_MSC_clusters <- quickCluster(TERT_MSC_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
TERT_7dOb_clusters <- quickCluster(TERT_7dOb_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
TERT_7dAd_clusters <- quickCluster(TERT_7dAd_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
TERT_7dOb_4dAd_clusters <- quickCluster(TERT_7dOb_4dAd_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
TERT_7dAd_4dOb_clusters <- quickCluster(TERT_7dAd_4dOb_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())

# Compute scaling factors
TERT_MSC_sce <- computeSumFactors(TERT_MSC_sce, min.mean=0.1, cluster=TERT_MSC_clusters)
TERT_7dOb_sce <- computeSumFactors(TERT_7dOb_sce, min.mean=0.1, cluster=TERT_7dOb_clusters)
TERT_7dAd_sce <- computeSumFactors(TERT_7dAd_sce, min.mean=0.1, cluster=TERT_7dAd_clusters)
TERT_7dOb_4dAd_sce <- computeSumFactors(TERT_7dOb_4dAd_sce, min.mean=0.1, cluster=TERT_7dOb_4dAd_clusters)
TERT_7dAd_4dOb_sce <- computeSumFactors(TERT_7dAd_4dOb_sce, min.mean=0.1, cluster=TERT_7dAd_4dOb_clusters)

# Normalize the counts
TERT_MSC_sce <- logNormCounts(TERT_MSC_sce)
TERT_7dOb_sce <- logNormCounts(TERT_7dOb_sce)
TERT_7dAd_sce <- logNormCounts(TERT_7dAd_sce)
TERT_7dOb_4dAd_sce <- logNormCounts(TERT_7dOb_4dAd_sce)
TERT_7dAd_4dOb_sce <- logNormCounts(TERT_7dAd_4dOb_sce)

## Calculate doublet scores. NOTE: THIS STEP IS NON-DETERMINISTI - RESULTS VARY FROM RUN TO RUN
TERT_MSC_sce$DoubletScore <- scDblFinder::computeDoubletDensity(TERT_MSC_sce, BSPARAM=IrlbaParam())
TERT_7dOb_sce$DoubletScore <- scDblFinder::computeDoubletDensity(TERT_7dOb_sce, BSPARAM=IrlbaParam())
TERT_7dAd_sce$DoubletScore <- scDblFinder::computeDoubletDensity(TERT_7dAd_sce, BSPARAM=IrlbaParam())
TERT_7dOb_4dAd_sce$DoubletScore <- scDblFinder::computeDoubletDensity(TERT_7dOb_4dAd_sce, BSPARAM=IrlbaParam())
TERT_7dAd_4dOb_sce$DoubletScore <- scDblFinder::computeDoubletDensity(TERT_7dAd_4dOb_sce, BSPARAM=IrlbaParam())

### QC by clustering - Individual datasets. NOTE: RERUNNING THE CODE MAY PRODUCE DIFFERENT CLUSTERS AND CLUSTER LABELS DUE TO NON-DETERMINISTIC STEPS
### RATIONaLE: Low quality nuclei may be included due to threshold effects. Deep clustering can help to reveal if there are groups of low quality nuclei that does not mix with the remaining nuclei, and thus can be removed.
## Create Seurat objects
TERT_MSC_seurat <- as.Seurat(TERT_MSC_sce, counts = "counts", data = "logcounts")
TERT_7dOb_seurat <- as.Seurat(TERT_7dOb_sce, counts = "counts", data = "logcounts")
TERT_7dAd_seurat <- as.Seurat(TERT_7dAd_sce, counts = "counts", data = "logcounts")
TERT_7dOb_4dAd_seurat <- as.Seurat(TERT_7dOb_4dAd_sce, counts = "counts", data = "logcounts")
TERT_7dAd_4dOb_seurat <- as.Seurat(TERT_7dAd_4dOb_sce, counts = "counts", data = "logcounts")

## TERT_MSC
# Iteration 1 - Use all genes
VariableFeatures(TERT_MSC_seurat) <- rownames(TERT_MSC_seurat)
TERT_MSC_seurat <- ScaleData(TERT_MSC_seurat)
TERT_MSC_seurat <- RunPCA(TERT_MSC_seurat)
TERT_MSC_seurat <- RunUMAP(TERT_MSC_seurat, dims=1:20, reduction="pca")
TERT_MSC_seurat <- FindNeighbors(object = TERT_MSC_seurat, dims = 1:20, reduction = "pca")
TERT_MSC_seurat <- FindClusters(TERT_MSC_seurat, resolution = 6, algorithm = 1)
VlnPlot(TERT_MSC_seurat, c("sum", "detected", "DoubletScore","subsets_Mito_percent"), pt.size=0)
# remove cell clusters with high doublet scores, high mitochondrial percents or low counts/genes
TERT_MSC_seurat <- subset(TERT_MSC_seurat, cells=rownames(TERT_MSC_seurat@meta.data[ !(TERT_MSC_seurat@meta.data$seurat_clusters %in% c(12,17,25,27,31)),]))

## TERT_7dOb
# Iteration 1 - Use all genes
VariableFeatures(TERT_7dOb_seurat) <- rownames(TERT_7dOb_seurat)
TERT_7dOb_seurat <- ScaleData(TERT_7dOb_seurat)
TERT_7dOb_seurat <- RunPCA(TERT_7dOb_seurat)
TERT_7dOb_seurat <- RunUMAP(TERT_7dOb_seurat, dims=1:20, reduction="pca")
TERT_7dOb_seurat <- FindNeighbors(object = TERT_7dOb_seurat, dims = 1:20, reduction = "pca")
TERT_7dOb_seurat <- FindClusters(TERT_7dOb_seurat, resolution = 6, algorithm = 1)
VlnPlot(TERT_7dOb_seurat, c("sum", "detected", "DoubletScore","subsets_Mito_percent"), pt.size=0)
# remove cell clusters with high doublet scores, high mitochondrial percents or low counts/genes
#NONE

## TERT_7dAd
# Iteration 1 - Use all genes
VariableFeatures(TERT_7dAd_seurat) <- rownames(TERT_7dAd_seurat)
TERT_7dAd_seurat <- ScaleData(TERT_7dAd_seurat)
TERT_7dAd_seurat <- RunPCA(TERT_7dAd_seurat)
TERT_7dAd_seurat <- RunUMAP(TERT_7dAd_seurat, dims=1:20, reduction="pca")
TERT_7dAd_seurat <- FindNeighbors(object = TERT_7dAd_seurat, dims = 1:20, reduction = "pca")
TERT_7dAd_seurat <- FindClusters(TERT_7dAd_seurat, resolution = 6, algorithm = 1)
VlnPlot(TERT_7dAd_seurat, c("sum", "detected", "DoubletScore","subsets_Mito_percent"), pt.size=0)
# remove cell clusters with high doublet scores, high mitochondrial percents or low counts/genes
TERT_7dAd_seurat <- subset(TERT_7dAd_seurat, cells=rownames(TERT_7dAd_seurat@meta.data[ !(TERT_7dAd_seurat@meta.data$seurat_clusters %in% c(45)),]))

## TERT_7dOb_4dAd
# Iteration 1 - Use all genes
VariableFeatures(TERT_7dOb_4dAd_seurat) <- rownames(TERT_7dOb_4dAd_seurat)
TERT_7dOb_4dAd_seurat <- ScaleData(TERT_7dOb_4dAd_seurat)
TERT_7dOb_4dAd_seurat <- RunPCA(TERT_7dOb_4dAd_seurat)
TERT_7dOb_4dAd_seurat <- RunUMAP(TERT_7dOb_4dAd_seurat, dims=1:20, reduction="pca")
TERT_7dOb_4dAd_seurat <- FindNeighbors(object = TERT_7dOb_4dAd_seurat, dims = 1:20, reduction = "pca")
TERT_7dOb_4dAd_seurat <- FindClusters(TERT_7dOb_4dAd_seurat, resolution = 6, algorithm = 1)
VlnPlot(TERT_7dOb_4dAd_seurat, c("sum", "detected", "DoubletScore","subsets_Mito_percent"), pt.size=0)
# remove cell clusters with high doublet scores, high mitochondrial percents or low counts/genes
TERT_7dOb_4dAd_seurat <- subset(TERT_7dOb_4dAd_seurat, cells=rownames(TERT_7dOb_4dAd_seurat@meta.data[ !(TERT_7dOb_4dAd_seurat@meta.data$seurat_clusters %in% c(0, 26, 53)),]))

## TERT_7dAd_4dOb
# Iteration 1 - Use all genes
VariableFeatures(TERT_7dAd_4dOb_seurat) <- rownames(TERT_7dAd_4dOb_seurat)
TERT_7dAd_4dOb_seurat <- ScaleData(TERT_7dAd_4dOb_seurat)
TERT_7dAd_4dOb_seurat <- RunPCA(TERT_7dAd_4dOb_seurat)
TERT_7dAd_4dOb_seurat <- RunUMAP(TERT_7dAd_4dOb_seurat, dims=1:20, reduction="pca")
TERT_7dAd_4dOb_seurat <- FindNeighbors(object = TERT_7dAd_4dOb_seurat, dims = 1:20, reduction = "pca")
TERT_7dAd_4dOb_seurat <- FindClusters(TERT_7dAd_4dOb_seurat, resolution = 6, algorithm = 1)
VlnPlot(TERT_7dAd_4dOb_seurat, c("sum", "detected", "DoubletScore","subsets_Mito_percent"), pt.size=0)
# remove cell clusters with high doublet scores, high mitochondrial percents or low counts/genes
# NONE

### QC by clustering - Combined NOTE: RERUNNING THE CODE MAY PRODUCE DIFFERENT CLUSTERS AND CLUSTER LABELS DUE TO NON-DETERMINISTIC STEPS
# Merge the datasets
sc_seurat <- merge(TERT_MSC_seurat, y=TERT_7dOb_seurat)
sc_seurat <- merge(sc_seurat, y=TERT_7dAd_seurat)
sc_seurat <- merge(sc_seurat, y=TERT_7dOb_4dAd_seurat)
sc_seurat <- merge(sc_seurat, y=TERT_7dAd_4dOb_seurat)

sc_seurat$Dataset <- 'none'
sc_seurat@meta.data[grep("TERT_MSC",rownames(sc_seurat@meta.data)),'Dataset'] <- 'TERT_MSC'
sc_seurat@meta.data[grep("TERT_7dOb",rownames(sc_seurat@meta.data)),'Dataset'] <- 'TERT_7dOb'
sc_seurat@meta.data[grep("TERT_7dAd",rownames(sc_seurat@meta.data)),'Dataset'] <- 'TERT_7dAd'
sc_seurat@meta.data[grep("TERT_7dOb_4dAd",rownames(sc_seurat@meta.data)),'Dataset'] <- 'TERT_7dOb_4dAd'
sc_seurat@meta.data[grep("TERT_7dAd_4dOb",rownames(sc_seurat@meta.data)),'Dataset'] <- 'TERT_7dAd_4dOb'

# Iteration 1. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
VariableFeatures(sc_seurat) <- rownames(sc_seurat)
sc_seurat <- ScaleData(sc_seurat)
sc_seurat <- RunPCA(sc_seurat)
sc_seurat <- RunUMAP(sc_seurat, dims=1:20, reduction="pca")
sc_seurat <- FindNeighbors(object = sc_seurat, dims = 1:20, reduction = "pca")
sc_seurat <- FindClusters(sc_seurat, resolution = 5, algorithm = 1)
VlnPlot(sc_seurat, c("sum", "detected", "DoubletScore","subsets_Mito_percent"), pt.size=0)
#NONE

#Subset each SCE object to the identified high-quality droplets. 
TERT_MSC_sce <- TERT_MSC_sce[,which(colnames(TERT_MSC_sce) %in% colnames(sc_seurat))]
TERT_7dOb_sce <- TERT_7dOb_sce[,which(colnames(TERT_7dOb_sce) %in% colnames(sc_seurat))]
TERT_7dAd_sce <- TERT_7dAd_sce[,which(colnames(TERT_7dAd_sce) %in% colnames(sc_seurat))]
TERT_7dOb_4dAd_sce <- TERT_7dOb_4dAd_sce[,which(colnames(TERT_7dOb_4dAd_sce) %in% colnames(sc_seurat))]
TERT_7dAd_4dOb_sce <- TERT_7dAd_4dOb_sce[,which(colnames(TERT_7dAd_4dOb_sce) %in% colnames(sc_seurat))]

#Normalize the datasets. NOTE: THIS STEP IS NON-DETERMINISTIC - RESULTS VARY FROM RUN TO RUN
# Cluster the nuclei in each dataset
TERT_MSC_clusters <- quickCluster(TERT_MSC_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
TERT_7dOb_clusters <- quickCluster(TERT_7dOb_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
TERT_7dAd_clusters <- quickCluster(TERT_7dAd_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
TERT_7dOb_4dAd_clusters <- quickCluster(TERT_7dOb_4dAd_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
TERT_7dAd_4dOb_clusters <- quickCluster(TERT_7dAd_4dOb_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())

# Compute scaling factors
TERT_MSC_sce <- computeSumFactors(TERT_MSC_sce, min.mean=0.1, cluster=TERT_MSC_clusters)
TERT_7dOb_sce <- computeSumFactors(TERT_7dOb_sce, min.mean=0.1, cluster=TERT_7dOb_clusters)
TERT_7dAd_sce <- computeSumFactors(TERT_7dAd_sce, min.mean=0.1, cluster=TERT_7dAd_clusters)
TERT_7dOb_4dAd_sce <- computeSumFactors(TERT_7dOb_4dAd_sce, min.mean=0.1, cluster=TERT_7dOb_4dAd_clusters)
TERT_7dAd_4dOb_sce <- computeSumFactors(TERT_7dAd_4dOb_sce, min.mean=0.1, cluster=TERT_7dAd_4dOb_clusters)

# Normalize the counts
TERT_MSC_sce <- logNormCounts(TERT_MSC_sce)
TERT_7dOb_sce <- logNormCounts(TERT_7dOb_sce)
TERT_7dAd_sce <- logNormCounts(TERT_7dAd_sce)
TERT_7dOb_4dAd_sce <- logNormCounts(TERT_7dOb_4dAd_sce)
TERT_7dAd_4dOb_sce <- logNormCounts(TERT_7dAd_4dOb_sce)

#Reduce batch effects by rescaling across datasets
# Normalize across samples
rescaled <- batchelor::multiBatchNorm(
  TERT_MSC_sce, 
  TERT_7dOb_sce,
  TERT_7dAd_sce, 
  TERT_7dOb_4dAd_sce,
  TERT_7dAd_4dOb_sce
)

### Merge all the data and embed
# Create seurat objects
TERT_MSC_seurat <- as.Seurat(rescaled[[1]], counts = "counts", data = "logcounts")
TERT_7dOb_seurat <- as.Seurat(rescaled[[2]], counts = "counts", data = "logcounts")
TERT_7dAd_seurat <- as.Seurat(rescaled[[3]], counts = "counts", data = "logcounts")
TERT_7dOb_4dAd_seurat <- as.Seurat(rescaled[[4]], counts = "counts", data = "logcounts")
TERT_7dAd_4dOb_seurat <- as.Seurat(rescaled[[5]], counts = "counts", data = "logcounts")

# Merge the datasets
sc_seurat <- merge(TERT_MSC_seurat, y=TERT_7dOb_seurat)
sc_seurat <- merge(sc_seurat, y=TERT_7dAd_seurat)
sc_seurat <- merge(sc_seurat, y=TERT_7dOb_4dAd_seurat)
sc_seurat <- merge(sc_seurat, y=TERT_7dAd_4dOb_seurat)

sc_seurat$Dataset <- 'none'
sc_seurat@meta.data[grep("TERT_MSC",rownames(sc_seurat@meta.data)),'Dataset'] <- 'TERT_MSC'
sc_seurat@meta.data[grep("TERT_7dOb",rownames(sc_seurat@meta.data)),'Dataset'] <- 'TERT_7dOb'
sc_seurat@meta.data[grep("TERT_7dAd",rownames(sc_seurat@meta.data)),'Dataset'] <- 'TERT_7dAd'
sc_seurat@meta.data[grep("TERT_7dOb_4dAd",rownames(sc_seurat@meta.data)),'Dataset'] <- 'TERT_7dOb_4dAd'
sc_seurat@meta.data[grep("TERT_7dAd_4dOb",rownames(sc_seurat@meta.data)),'Dataset'] <- 'TERT_7dAd_4dOb'

# Make clusters on sc_seurat 
sc_seurat <- ScaleData(sc_seurat)
sc_seurat <- FindVariableFeatures(sc_seurat)
sc_seurat <- RunPCA(sc_seurat)
ElbowPlot(sc_seurat)
sc_seurat <- RunUMAP(sc_seurat, dims=1:19, reduction="pca") 
sc_seurat <- FindNeighbors(object = sc_seurat, dims = 1:14, reduction = "pca")
sc_seurat <- FindClusters(sc_seurat, resolution = 0.04, algorithm = 3)
DimPlot(sc_seurat, group.by="Dataset")

# Incorporate Ad and Ob scaled / percentage data
Matrix <- data.frame(as.matrix(GetAssayData(sc_seurat,slot="data")))
Matrix_scale <- data.frame(as.matrix(GetAssayData(sc_seurat,slot="scale.data")))

# gene list of Ob and Ad specific gene-markers
Markers <- readRDS("sc_gene_markers_1.Rds")
colnames(Markers)[1] <- "Symbol"

ScaleMarker <- data.frame(
  'Ob_scaled' = colSums(Matrix_scale[rownames(Matrix_scale) %in% Markers[Markers$Ob ==1,'Symbol',],]),
  'Ad_scaled' = colSums(Matrix_scale[rownames(Matrix_scale) %in% Markers[Markers$Ad ==1,'Symbol',],]))

sc_seurat$Ob_scaled <- ScaleMarker$Ob_scaled
sc_seurat$Ad_scaled <- ScaleMarker$Ad_scaled
#saveRDS(sc_seurat, "sc_seurat_diff_scales.Rds")

  ### snRNA-seq seurat object
Dataset <- as.data.frame(as.matrix(Matrix::readMM("Parse_transdiff.mtx")))

genes <- read.csv("all_genes_parse.csv") 
cells <- read.csv("cell_metadata_parse.csv") 

colnames(Dataset) <- genes$gene_name

cells$sublibrary <- NA
cells[1:dim(cells[grep("s1",cells$bc_wells),])[1],]$sublibrary <- "lib_1"
cells[((dim(cells[grep("s1",cells$bc_wells),])[1])+1):ncol(Dataset),]$sublibrary <- "lib_2"
cells <- na.omit(cells)

cells$ID <- NA
cells$ID <- paste(cells$bc1_well, cells$bc2_well, cells$bc3_well, cells$sublibrary, sep="_")

Dataset <- Dataset[1:nrow(cells),]
rownames(Dataset) <- cells$ID

Dataset <- t(Dataset)
Dataset <- as.data.frame(Dataset)

##add sample name
cells[cells$sample %in% "xcond_1",]$sample <- "pMSC_1"
cells[cells$sample %in% "xcond_2",]$sample <- "pMSC_2"
cells[cells$sample %in% "xcond_3",]$sample <- "StartOb"
cells[cells$sample %in% "xcond_4",]$sample <- "StartAd"

cells$ID2 <- paste(cells$sample, cells$ID, sep="_")

#Remove unnecessary samples
colnames(Dataset) <- cells$ID2
Dataset <- Dataset[,colnames(Dataset) %!in% c(colnames(Dataset[,grep("pMSC_1", colnames(Dataset))]), colnames(Dataset[,grep("pMSC_2", colnames(Dataset))]), colnames(Dataset[,grep("StartAd", colnames(Dataset))]))]

## Quality control
## Create SingleCellExperiment objects to use DropUtils, scran, scater, etc.
Dataset_sce <- SingleCellExperiment(assays = list(counts=Dataset))

## Calculate QC parameters (throws a warning, that can be ignored)
Dataset_sce <- addPerCellQC(Dataset_sce, subsets=list(Mito=grep(pattern = "^MT-", x = rownames(Dataset_sce), value = FALSE)))

## Threshold filtering of droplets in each dataset
## REMOVE: Droplets with more than 20% mitochondrial reads, less than 1000 UMIs, less than 500 genes or extremely high ratio between counts and genes (low complexity))
Dataset_sce <- Dataset_sce[,!(Dataset_sce$subsets_Mito_percent > 20)] 
Dataset_sce <- Dataset_sce[,!(Dataset_sce$sum/Dataset_sce$detected > 15)] 
Dataset_sce <- Dataset_sce[,(Dataset_sce$sum >= 100 & Dataset_sce$detected >= 100)] 
dim(Dataset_sce) 


## Automatic filtering droplets in each dataset using PCA across all QC metrics
# Calculate outliers
Dataset_sce <- logNormCounts(Dataset_sce)

Dataset_sce <- runColDataPCA(Dataset_sce,variables=list("sum","detected","subsets_Mito_percent"),outliers=T)

table(Dataset_sce$outlier)
# Remove outliers with PCA
Dataset_sce <- Dataset_sce[ ,!Dataset_sce$outlier] 

## Threshold filtering of genes in each dataset
## REMOVE: Genes expressed in less than 10 nuclei in all datasets
# Find lowly expressed genes and get the intersection
Dataset_low <- names(which(nexprs(Dataset_sce, byrow=T) <= 1))

Low <- Reduce(intersect, list(Dataset_low))

# Remove lowly expressed genes
Dataset_sce <- Dataset_sce[ which(!(rownames(Dataset_sce) %in% Low)),] 

## Filtering genes based on biotype and transcript level support (not done yet)
## REMOVE: Non-protein coding genes, keep the one with high transcript level support (tsl)
Genes <- read.delim("mart_export_human_Ens_GRCh38.txt")
Genes <- Genes[ Genes$Gene.name %in% rownames(P1_sce),]
Genes <- Genes[ grep("tsl1", Genes$Transcript.support.level..TSL.),]
Genes <- Genes[ Genes$Gene.name !=  "",]
Genes <- Genes[ Genes$Transcript.type == "protein_coding",]
Genes <- Genes[ Genes$Gene.type == "protein_coding",]
Genes <- Genes[ !is.na(Genes$Gene.name),]
Genes <- Genes[ duplicated(Genes$Gene.stable.ID) == F,] 

# Filter the sce objects and ambient genes
Dataset_sce <- Dataset_sce[ which(rownames(Dataset_sce) %in% Genes$Gene.name),] # 16094 genes left

## Normalize the count matrices #
# Cluster each data. NOTE: THIS STEP IS NON-DETERMINISTIC - RESULTS VARY FROM RUN TO RUN
Dataset_clusters <- quickCluster(Dataset_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())

# Compute scaling factors
Dataset_sce <- computeSumFactors(Dataset_sce, min.mean=0.1, cluster=Dataset_clusters)

# Normalize the counts
Dataset_sce <- logNormCounts(Dataset_sce)

## Calculate doublet scores. NOTE: THIS STEP IS NON-DETERMINISTI - RESULTS VARY FROM RUN TO RUN
Dataset_sce$DoubletScore <- scDblFinder::computeDoubletDensity(Dataset_sce, BSPARAM=IrlbaParam())

### QC by clustering - Individual datasets. NOTE: RERUNNING THE CODE MAY PRODUCE DIFFERENT CLUSTERS AND CLUSTER LABELS DUE TO NON-DETERMINISTIC STEPS
### RATIONaLE: Low quality nuclei may be included due to threshold effects. Deep clustering can help to reveal if there are groups of low quality nuclei that does not mix with the remaining nuclei, and thus can be removed.
## Create Seurat objects
Dataset_seurat <- as.Seurat(Dataset_sce, counts = "counts", data = "logcounts")

## P1
# Iteration 1 - Use all genes
VariableFeatures(Dataset_seurat) <- rownames(Dataset_seurat)
Dataset_seurat <- ScaleData(Dataset_seurat)
Dataset_seurat <- RunPCA(Dataset_seurat)
Dataset_seurat <- RunUMAP(Dataset_seurat, dims=1:20, reduction="pca")
Dataset_seurat <- FindNeighbors(object = Dataset_seurat, dims = 1:20, reduction = "pca")
Dataset_seurat <- FindClusters(Dataset_seurat, resolution = 6, algorithm = 1)
# remove cell clusters with high doublet scores, high mitochondrial percents or low counts/genes
Dataset_seurat <- subset(Dataset_seurat, cells=rownames(Dataset_seurat@meta.data[ !(Dataset_seurat@meta.data$seurat_clusters %in% c(5,48)),]))

Dataset_seurat$Dataset <- 'none'
Dataset_seurat@meta.data[grep("StartOb",rownames(Dataset_seurat@meta.data)),'Dataset'] <- 'StartOb'

Dataset_sce <- Dataset_sce[,which(colnames(Dataset_sce) %in% colnames(Dataset_seurat))]
Dataset_clusters <- quickCluster(Dataset_sce, use.ranks=FALSE, BSPARAM=IrlbaParam()) 
Dataset_sce <- computeSumFactors(Dataset_sce, min.mean=0.1, cluster=Dataset_clusters)
Dataset_sce <- logNormCounts(Dataset_sce)

Dataset_seurat <- as.Seurat(Dataset_sce, counts = "counts", data = "logcounts")

Dataset_seurat$Dataset <- 'none'
Dataset_seurat@meta.data[grep("StartOb",rownames(Dataset_seurat@meta.data)),'Dataset'] <- 'StartOb'

Dataset_seurat <- ScaleData(Dataset_seurat, features=rownames(Dataset_seurat))
Dataset_seurat <- FindVariableFeatures(Dataset_seurat) 
Dataset_seurat <- RunPCA(Dataset_seurat)
Dataset_seurat <- RunUMAP(Dataset_seurat, dims=c(1:20), reduction="pca")
Dataset_seurat <- FindNeighbors(object = Dataset_seurat, dims = 1:2, reduction = "umap")
Dataset_seurat <- FindClusters(Dataset_seurat, resolution = 0.05, algorithm = 3)
#saveRDS(Dataset_seurat, "Hybrid_seurat_3.Rds")



















