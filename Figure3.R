### Figure 3
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

# load the following objects (see DataPreparation.R for more information on how these were generated)
seurat_all <- readRDS("sc_seurat_diff_scales.Rds")

### Figure 3B
# UMAP-plot of seurat object
DimPlot(seurat_all, group.by="Dataset")

### Figure 3C
# Violin plot of Osteoblast (ALPL & LMCD1), Adipocyte (PLIN1 & ADIPOQ), and stem cell (ENG, THY1, and LEPR) markers within each dataset
seurat_all$Dataset <- factor(seurat_all$Dataset, levels = c("TERT_7dOb_4dAd","TERT_7dOb","TERT_MSC","TERT_7dAd","TERT_7dAd_4dOb"))
VlnPlot(seurat_all, c("ALPL", "LMCD1"), group.by="Dataset")
VlnPlot(seurat_all, c("PLIN1","ADIPOQ"), group.by="Dataset")
VlnPlot(seurat_all, c("THY1","LEPR"), group.by="Dataset")
VlnPlot(seurat_all, c("ENG"), group.by="Dataset")

### Figure 3D
library(Seurat)
FeaturePlot(seurat_all, "Ad_scaled")
FeaturePlot(seurat_all, "Ob_scaled")

### Figure 3E
# calculate pseudotime values for each dataset using monocle2

#7dAD
seurat_7dAd <- subset(seurat_all, cells=rownames(seurat_all@meta.data[seurat_all@meta.data$Dataset %in% "TERT_7dAd",]))

Tmp <- data.frame(gene_short_name = rownames(seurat_7dAd))
rownames(Tmp) <- rownames(seurat_7dAd)
fd <- new("AnnotatedDataFrame", data = Tmp)
pd <- new("AnnotatedDataFrame", data = seurat_7dAd@meta.data)
Mncl.std <- newCellDataSet(GetAssayData(seurat_7dAd, slot = "scale.data"),phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = gaussianff()) #Call uninormal() instead of gaussianff()

# Use all variable genes
genes <- VariableFeatures(seurat_7dAd)
Mncl.std  <- setOrderingFilter(Mncl.std , genes)

#Order cells
Mncl.std  <- reduceDimension(Mncl.std , max_components = 2, method = 'DDRTree', norm_method = "none", verbose = T, pseudo_expr = 0)
Mncl.std  <- orderCells(Mncl.std )

#Plot Mncl
plot_cell_trajectory(Mncl.std , color_by = "State")

#Choose root state
Mncl.std  <- orderCells(Mncl.std , root_state = 1) #choose root states

#Extract pseudotime information
Tmp <- as.data.frame(cbind(Mncl.std$State,Mncl.std$Dataset,Mncl.std$Pseudotime))
Tmp$Cell <- colnames(Mncl.std)
colnames(Tmp) <- c("State","Dataset","Pseudotime","Cell")
Tmp2 <- as.data.frame(seurat_7dAd@meta.data)
Tmp2$Cell <- rownames(Tmp2)
Tmp2 <- merge(Tmp2, Tmp[,c("Pseudotime","Cell","State")], sort=F, by="Cell")

#Add information to seurat
seurat_7dAd$Pseudotime <- as.numeric(Tmp2$Pseudotime)
seurat_7dAd$State <- as.numeric(Tmp2$State)
# 1st final plot for Figure 3E
FeaturePlot(seurat_7dAd, "Pseudotime")

#7dOb
seurat_7dOb <- subset(seurat_all, cells=rownames(seurat_all@meta.data[seurat_all@meta.data$Dataset %in% "TERT_7dOb",]))

Tmp <- data.frame(gene_short_name = rownames(seurat_7dOb))
rownames(Tmp) <- rownames(seurat_7dOb)
fd <- new("AnnotatedDataFrame", data = Tmp)
pd <- new("AnnotatedDataFrame", data = seurat_7dOb@meta.data)
Mncl.std <- newCellDataSet(GetAssayData(seurat_7dOb, slot = "scale.data"),phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = gaussianff()) #Call uninormal() instead of gaussianff()

# Use all variable genes
genes <- VariableFeatures(seurat_7dOb)
Mncl.std  <- setOrderingFilter(Mncl.std , genes)

#Order cells
Mncl.std  <- reduceDimension(Mncl.std , max_components = 2, method = 'DDRTree', norm_method = "none", verbose = T, pseudo_expr = 0)
Mncl.std  <- orderCells(Mncl.std )

#Plot Mncl
plot_cell_trajectory(Mncl.std , color_by = "State")
plot_cell_trajectory(Mncl.std , color_by = "Pseudotime")

#Choose root state
Mncl.std  <- orderCells(Mncl.std , root_state = 4) #choose root states
#Reverse the order
pData(Mncl.std)$Pseudotime <- max(pData(Mncl.std)$Pseudotime) - pData(Mncl.std)$Pseudotime # to reverse

#Extract pseudotime information
Tmp <- as.data.frame(cbind(Mncl.std$State,Mncl.std$Dataset,Mncl.std$Pseudotime))
Tmp$Cell <- colnames(Mncl.std)
colnames(Tmp) <- c("State","Dataset","Pseudotime","Cell")
Tmp2 <- as.data.frame(seurat_7dOb@meta.data)
Tmp2$Cell <- rownames(Tmp2)
Tmp2 <- merge(Tmp2, Tmp[,c("Pseudotime","Cell","State")], sort=F, by="Cell")

#Add information to seurat
seurat_7dOb$Pseudotime <- as.numeric(Tmp2$Pseudotime)
seurat_7dOb$State <- as.numeric(Tmp2$State)
# 2nd final plot for Figure 3E
FeaturePlot(seurat_7dOb, seurat_7dOb, "Pseudotime")

#7dOb_4dAd
seurat_7dOb_4dAd <- subset(seurat_all, cells=rownames(seurat_all@meta.data[seurat_all@meta.data$Dataset %in% "TERT_7dOb_4dAd",]))

Tmp <- data.frame(gene_short_name = rownames(seurat_7dOb_4dAd))
rownames(Tmp) <- rownames(seurat_7dOb_4dAd)
fd <- new("AnnotatedDataFrame", data = Tmp)
pd <- new("AnnotatedDataFrame", data = seurat_7dOb_4dAd@meta.data)
Mncl.std <- newCellDataSet(GetAssayData(seurat_7dOb_4dAd, slot = "scale.data"),phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = gaussianff()) #Call uninormal() instead of gaussianff()

# Use all variable genes
genes <- VariableFeatures(seurat_7dOb_4dAd)
Mncl.std  <- setOrderingFilter(Mncl.std , genes)

#Order cells
Mncl.std  <- reduceDimension(Mncl.std , max_components = 2, method = 'DDRTree', norm_method = "none", verbose = T, pseudo_expr = 0)
Mncl.std  <- orderCells(Mncl.std )

#Plot Mncl
plot_cell_trajectory(Mncl.std , color_by = "State")

#Choose root state
Mncl.std  <- orderCells(Mncl.std , root_state = 1) #choose root states
#Reverse the order
pData(Mncl.std)$Pseudotime <- max(pData(Mncl.std)$Pseudotime) - pData(Mncl.std)$Pseudotime # to reverse

#Extract pseudotime information
Tmp <- as.data.frame(cbind(Mncl.std$State,Mncl.std$Dataset,Mncl.std$Pseudotime))
Tmp$Cell <- colnames(Mncl.std)
colnames(Tmp) <- c("State","Dataset","Pseudotime","Cell")
Tmp2 <- as.data.frame(seurat_7dOb_4dAd@meta.data)
Tmp2$Cell <- rownames(Tmp2)
Tmp2 <- merge(Tmp2, Tmp[,c("Pseudotime","Cell","State")], sort=F, by="Cell")

#Add information to seurat
seurat_7dOb_4dAd$Pseudotime <- as.numeric(Tmp2$Pseudotime)
seurat_7dOb_4dAd$State <- as.numeric(Tmp2$State)
# 3rd final plot for Figure 3E
FeaturePlot(seurat_7dOb_4dAd, "Pseudotime")

#7dAd_4dOb
seurat_7dAd_4dOb <- subset(seurat_all, cells=rownames(seurat_all@meta.data[seurat_all@meta.data$Dataset %in% "TERT_7dAd_4dOb",]))

Tmp <- data.frame(gene_short_name = rownames(seurat_7dAd_4dOb))
rownames(Tmp) <- rownames(seurat_7dAd_4dOb)
fd <- new("AnnotatedDataFrame", data = Tmp)
pd <- new("AnnotatedDataFrame", data = seurat_7dAd_4dOb@meta.data)
Mncl.std <- newCellDataSet(GetAssayData(seurat_7dAd_4dOb, slot = "scale.data"),phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = gaussianff()) #Call uninormal() instead of gaussianff()

# Use all variable genes
genes <- VariableFeatures(seurat_7dAd_4dOb)
Mncl.std  <- setOrderingFilter(Mncl.std , genes)

#Order cells
Mncl.std  <- reduceDimension(Mncl.std , max_components = 2, method = 'DDRTree', norm_method = "none", verbose = T, pseudo_expr = 0)
Mncl.std  <- orderCells(Mncl.std )

#Plot Mncl
plot_cell_trajectory(Mncl.std , color_by = "State")

#Choose root state
Mncl.std  <- orderCells(Mncl.std , root_state = 1) #choose root states
#Reverse the order
pData(Mncl.std)$Pseudotime <- max(pData(Mncl.std)$Pseudotime) - pData(Mncl.std)$Pseudotime # to reverse

#Extract pseudotime information
Tmp <- as.data.frame(cbind(Mncl.std$State,Mncl.std$Dataset,Mncl.std$Pseudotime))
Tmp$Cell <- colnames(Mncl.std)
colnames(Tmp) <- c("State","Dataset","Pseudotime","Cell")
Tmp2 <- as.data.frame(seurat_7dAd_4dOb@meta.data)
Tmp2$Cell <- rownames(Tmp2)
Tmp2 <- merge(Tmp2, Tmp[,c("Pseudotime","Cell","State")], sort=F, by="Cell")

#Add information to seurat
seurat_7dAd_4dOb$Pseudotime <- as.numeric(Tmp2$Pseudotime)
seurat_7dAd_4dOb$State <- as.numeric(Tmp2$State)
# 4th final plot for Figure 3E
FeaturePlot(seurat_7dAd_4dOb, "Pseudotime")

## Find expressed genes (in at least 10% of the cells) across direct and their trans-differentiated datasets, towards Ob
sca <-  FromMatrix(as.matrix(seurat_7dOb@assays$RNA@data), seurat_7dOb@meta.data)
expressed_genes1 <- names(which(freq(sca) >= 0.1))
sca <-  FromMatrix(as.matrix(seurat_4dOb@assays$RNA@data), seurat_4dOb@meta.data)
expressed_genes2 <- names(which(freq(sca) >= 0.1))
expressed_genes <- unique(c(expressed_genes1, expressed_genes2))

# 7dOb
Pseudotime <- data.frame(curve1 = seurat_7dOb$Pseudotime)
Pseudotime <- na.omit(Pseudotime)
Weights <- data.frame(TERT_7dOb = rep(1, nrow(Pseudotime)))
rownames(Weights) = rownames(Pseudotime)
Counts <- seurat_7dOb@assays$RNA@counts
Counts <- Counts[,colnames(Counts) %in% rownames(Pseudotime)]
Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
#Fit the model
fitGAM_7dOb <- fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 10, verbose = TRUE)

# 4dOb
Pseudotime <- data.frame(curve1 = seurat_4dOb$Pseudotime)
Pseudotime <- na.omit(Pseudotime)
Weights <- data.frame(TERT_7dOb_4dAd = rep(1, nrow(Pseudotime)))
rownames(Weights) = rownames(Pseudotime)
Counts <- seurat_4dOb@assays$RNA@counts
Counts <- Counts[,colnames(Counts) %in% rownames(Pseudotime)]
Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
#Fit the model
fitGAM_4dOb <- fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 10, verbose = TRUE)

## Find expressed genes (in at least 10% of the cells) across direct and their trans-differentiated datasets, towards Ad
sca <-  FromMatrix(as.matrix(seurat_7dAd@assays$RNA@data), seurat_7dAd@meta.data)
expressed_genes1 <- names(which(freq(sca) >= 0.1))
sca <-  FromMatrix(as.matrix(seurat_4dAd@assays$RNA@data), seurat_4dAd@meta.data)
expressed_genes2 <- names(which(freq(sca) >= 0.1))
expressed_genes <- unique(c(expressed_genes1, expressed_genes2))

# 7dAd
Pseudotime <- data.frame(curve1 = seurat_7dAd$Pseudotime)
Pseudotime <- na.omit(Pseudotime)
Weights <- data.frame(TERT_7dAd = rep(1, nrow(Pseudotime)))
rownames(Weights) = rownames(Pseudotime)
Counts <- seurat_7dAd@assays$RNA@counts
Counts <- Counts[,colnames(Counts) %in% rownames(Pseudotime)]
Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
#Fit the model
fitGAM_7dAd <- fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 10, verbose = TRUE)

# 4dAd
Pseudotime <- data.frame(curve1 = seurat_4dAd$Pseudotime)
Pseudotime <- na.omit(Pseudotime)
Weights <- data.frame(TERT_7dOb_4dAd = rep(1, nrow(Pseudotime)))
rownames(Weights) = rownames(Pseudotime)
Counts <- seurat_4dAd@assays$RNA@counts
Counts <- Counts[,colnames(Counts) %in% rownames(Pseudotime)]
Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
#Fit the model
fitGAM_4dAd <- fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 10, verbose = TRUE)

#save fitGAM objects
saveRDS(fitGAM_7dAd, "x7dAd_only_woProf_fitGAM_2")
saveRDS(fitGAM_4dAd, "x7dOb4dAd_only_woProf_fitGAM_2")
saveRDS(fitGAM_7dOb, "x7dOb_sce_fitGAM_new")
saveRDS(fitGAM_4dOb, "x7dAd4dOb_only_woProf_fitGAM_new")

### Figure 3F
# Hexbinplot of Ad scaled vs. Ob scaled for each dataset
#load object
seurat_scaled <- readRDS("sc_seurat_diff_scales.Rds")

#TERT_7dOb
x <- subset(seurat_scaled, cells=rownames(seurat_scaled@meta.data[seurat_scaled@meta.data$Dataset %in% "TERT_7dOb",]))
df <- as.data.frame(cbind(x$Ad_scaled, x$Ob_scaled))
colnames(df) <- c("Ad_scaled", "Ob_scaled")
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
x <- 70
hexbinplot(Ad_scaled ~ Ob_scaled,data=df, xlab="Ob scaled", ylab="Ad scaled", colramp=rf, xbins=x, aspect=1, ylim=c(-60,480), xlim=c(-70,150))

#TERT_7dOb_4dAd
x <- subset(seurat_scaled, cells=rownames(seurat_scaled@meta.data[seurat_scaled@meta.data$Dataset %in% "TERT_7dOb_4dAd",]))
df <- as.data.frame(cbind(x$Ad_scaled, x$Ob_scaled))
colnames(df) <- c("Ad_scaled", "Ob_scaled")
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
x <- 70
hexbinplot(Ad_scaled ~ Ob_scaled,data=df, xlab="Ob scaled", ylab="Ad scaled", colramp=rf, xbins=x, aspect=1, ylim=c(-60,480), xlim=c(-70,150))

#TERT_7dAd
x <- subset(seurat_scaled, cells=rownames(seurat_scaled@meta.data[seurat_scaled@meta.data$Dataset %in% "TERT_7dAd",]))
df <- as.data.frame(cbind(x$Ad_scaled, x$Ob_scaled))
colnames(df) <- c("Ad_scaled", "Ob_scaled")
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
x <- 45
hexbinplot(Ad_scaled ~ Ob_scaled,data=df, xlab="Ob scaled", ylab="Ad scaled", colramp=rf, xbins=x, aspect=1, ylim=c(-60,480), xlim=c(-70,150))

#TERT_7dAd_4dOb
x <- subset(seurat_scaled, cells=rownames(seurat_scaled@meta.data[seurat_scaled@meta.data$Dataset %in% "TERT_7dAd_4dOb",]))
df <- as.data.frame(cbind(x$Ad_scaled, x$Ob_scaled))
colnames(df) <- c("Ad_scaled", "Ob_scaled")
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
x <- 75
hexbinplot(Ad_scaled ~ Ob_scaled,data=df, xlab="Ob scaled", ylab="Ad scaled", colramp=rf, xbins=x, aspect=1, ylim=c(-60,480), xlim=c(-70,150))

### Figure 3G & 3J
## Smooth heatmap comparisons

# load objects
MSC_sce_7dAd <- readRDS("x7dAd_only_woProf_fitGAM_2")
MSC_sce_4dAd <- readRDS("x7dOb4dAd_only_woProf_fitGAM_2")
MSC_sce_7dOb <- readRDS("x7dOb_sce_fitGAM_new")
MSC_sce_4dOb <- readRDS("x7dAd4dOb_only_woProf_fitGAM_new")

Ass_7dAd <- as.data.frame(associationTest(MSC_sce_7dAd))
Ass_4dAd <- as.data.frame(associationTest(MSC_sce_4dAd))
Ass_7dOb <- as.data.frame(associationTest(MSC_sce_7dOb))
Ass_4dOb <- as.data.frame(associationTest(MSC_sce_4dOb))

Ass_7dAd[!complete.cases(Ass_7dAd$pvalue),'pvalue'] <- 1
Ass_4dAd[!complete.cases( Ass_4dAd$pvalue),'pvalue'] <- 1
Ass_7dOb[!complete.cases(Ass_7dOb$pvalue),'pvalue'] <- 1
Ass_4dOb[!complete.cases( Ass_4dOb$pvalue),'pvalue'] <- 1

#Insert p_value threshold
pvalue <- 0.001
Ass_7dAd_filter <- Ass_7dAd[Ass_7dAd$pvalue < pvalue,]
Ass_4dAd_filter <- Ass_4dAd[Ass_4dAd$pvalue < pvalue,]
Ass_7dOb_filter <- Ass_7dOb[Ass_7dOb$pvalue < pvalue,]
Ass_4dOb_filter <- Ass_4dOb[Ass_4dOb$pvalue < pvalue,]

Ass_7dAd_filter$genes <- rownames(Ass_7dAd_filter)
Ass_4dAd_filter$genes <- rownames(Ass_4dAd_filter)
Ass_7dOb_filter$genes <- rownames(Ass_7dOb_filter)
Ass_4dOb_filter$genes <- rownames(Ass_4dOb_filter)

#common genes across both datasets
filter_Ad <- merge(Ass_7dAd_filter, Ass_4dAd_filter, by="genes")$genes #2287
filter_Ob <- merge(Ass_7dOb_filter, Ass_4dOb_filter, by="genes")$genes #3621

# Make "not in" operator
`%!in%` <- Negate(`%in%`)

#Unique genes across datasets
#Ad
filter_7dAd_only <- Ass_7dAd_filter[Ass_7dAd_filter$genes %!in% filter_Ad,]$genes #685
filter_4dAd_only <- Ass_4dAd_filter[Ass_4dAd_filter$genes %!in% filter_Ad,]$genes #553

#Ob
filter_7dOb_only <- Ass_7dOb_filter[Ass_7dOb_filter$genes %!in% filter_Ob,]$genes #997
filter_4dOb_only <- Ass_4dOb_filter[Ass_4dOb_filter$genes %!in% filter_Ob,]$genes #1003

#Subset Ass_xx to only and common genes for each direction and library
Ass_7dAd_1 <- Ass_7dAd[rownames(Ass_7dAd) %in% filter_7dAd_only,]
Ass_7dAd_2 <- Ass_7dAd[rownames(Ass_7dAd) %in% filter_Ad,]
Ass_7dAd_3 <- Ass_7dAd[rownames(Ass_7dAd) %in% filter_4dAd_only,]

Ass_4dAd_1 <- Ass_4dAd[rownames(Ass_4dAd) %in% filter_7dAd_only,]
Ass_4dAd_2 <- Ass_4dAd[rownames(Ass_4dAd) %in% filter_Ad,]
Ass_4dAd_3 <- Ass_4dAd[rownames(Ass_4dAd) %in% filter_4dAd_only,]

Ass_7dOb_1 <- Ass_7dOb[rownames(Ass_7dOb) %in% filter_7dOb_only,]
Ass_7dOb_2 <- Ass_7dOb[rownames(Ass_7dOb) %in% filter_Ob,]
Ass_7dOb_3 <- Ass_7dOb[rownames(Ass_7dOb) %in% filter_4dOb_only,]

Ass_4dOb_1 <- Ass_4dOb[rownames(Ass_4dOb) %in% filter_7dOb_only,]
Ass_4dOb_2 <- Ass_4dOb[rownames(Ass_4dOb) %in% filter_Ob,]
Ass_4dOb_3 <- Ass_4dOb[rownames(Ass_4dOb) %in% filter_4dOb_only,]

#PredictSmooth
Smooth_MSC_7dAd_1 <- predictSmooth(MSC_sce_7dAd, gene = rownames(Ass_7dAd_1), tidy=F, n=100)
Smooth_MSC_7dAd_2 <- predictSmooth(MSC_sce_7dAd, gene = rownames(Ass_7dAd_2), tidy=F, n=100)
Smooth_MSC_7dAd_3 <- predictSmooth(MSC_sce_7dAd, gene = rownames(Ass_7dAd_3), tidy=F, n=100)

Smooth_MSC_4dAd_1 <- predictSmooth(MSC_sce_4dAd, gene = rownames(Ass_4dAd_1), tidy=F, n=100)
Smooth_MSC_4dAd_2 <- predictSmooth(MSC_sce_4dAd, gene = rownames(Ass_4dAd_2), tidy=F, n=100)
Smooth_MSC_4dAd_3 <- predictSmooth(MSC_sce_4dAd, gene = rownames(Ass_4dAd_3), tidy=F, n=100)

Smooth_MSC_7dOb_1 <- predictSmooth(MSC_sce_7dOb, gene = rownames(Ass_7dOb_1), tidy=F, n=100)
Smooth_MSC_7dOb_2 <- predictSmooth(MSC_sce_7dOb, gene = rownames(Ass_7dOb_2), tidy=F, n=100)
Smooth_MSC_7dOb_3 <- predictSmooth(MSC_sce_7dOb, gene = rownames(Ass_7dOb_3), tidy=F, n=100)

Smooth_MSC_4dOb_1 <- predictSmooth(MSC_sce_4dOb, gene = rownames(Ass_4dOb_1), tidy=F, n=100)
Smooth_MSC_4dOb_2 <- predictSmooth(MSC_sce_4dOb, gene = rownames(Ass_4dOb_2), tidy=F, n=100)
Smooth_MSC_4dOb_3 <- predictSmooth(MSC_sce_4dOb, gene = rownames(Ass_4dOb_3), tidy=F, n=100)

# Scaling
Smooth_MSC_7dAd_1 <- t(scale(t(Smooth_MSC_7dAd_1)))
Smooth_MSC_7dAd_2 <- t(scale(t(Smooth_MSC_7dAd_2)))
Smooth_MSC_7dAd_3 <- t(scale(t(Smooth_MSC_7dAd_3)))

Smooth_MSC_4dAd_1 <- t(scale(t(Smooth_MSC_4dAd_1)))
Smooth_MSC_4dAd_2 <- t(scale(t(Smooth_MSC_4dAd_2)))
Smooth_MSC_4dAd_3 <- t(scale(t(Smooth_MSC_4dAd_3)))

Smooth_MSC_7dOb_1 <- t(scale(t(Smooth_MSC_7dOb_1)))
Smooth_MSC_7dOb_2 <- t(scale(t(Smooth_MSC_7dOb_2)))
Smooth_MSC_7dOb_3 <- t(scale(t(Smooth_MSC_7dOb_3)))

Smooth_MSC_4dOb_1 <- t(scale(t(Smooth_MSC_4dOb_1)))
Smooth_MSC_4dOb_2 <- t(scale(t(Smooth_MSC_4dOb_2)))
Smooth_MSC_4dOb_3 <- t(scale(t(Smooth_MSC_4dOb_3)))

# Seriate the results based on the significanrt group
# filter_7dAd_only
Smooth_MSC_7dAd_1 <- Smooth_MSC_7dAd_1[ get_order(seriate(Smooth_MSC_7dAd_1, method="PCA_angle")),]
# filter_Ad - does not matter which data frame to chose 
Smooth_MSC_7dAd_2 <- Smooth_MSC_7dAd_2[ get_order(seriate(Smooth_MSC_7dAd_2, method="PCA_angle")),]
# filter_4dAd_only
Smooth_MSC_4dAd_3 <- Smooth_MSC_4dAd_3[ get_order(seriate(Smooth_MSC_4dAd_3, method="PCA_angle")),]

# filter_7dOb_only
Smooth_MSC_7dOb_1 <- Smooth_MSC_7dOb_1[ get_order(seriate(Smooth_MSC_7dOb_1, method="PCA_angle")),]
# filter_Ob - does not matter which data frame to chose 
Smooth_MSC_7dOb_2 <- Smooth_MSC_7dOb_2[ get_order(seriate(Smooth_MSC_7dOb_2, method="PCA_angle")),]
# filter_4dOb_only
Smooth_MSC_4dOb_3 <- Smooth_MSC_4dOb_3[ get_order(seriate(Smooth_MSC_4dOb_3, method="PCA_angle")),]

# Arrange rows in the same order
Smooth_MSC_4dAd_1 <- Smooth_MSC_4dAd_1[rownames(Smooth_MSC_7dAd_1),]
Smooth_MSC_4dAd_2 <- Smooth_MSC_4dAd_2[rownames(Smooth_MSC_7dAd_2),]
Smooth_MSC_7dAd_3 <- Smooth_MSC_7dAd_3[rownames(Smooth_MSC_4dAd_3),]

Smooth_MSC_4dOb_1 <- Smooth_MSC_4dOb_1[rownames(Smooth_MSC_7dOb_1),]
Smooth_MSC_4dOb_2 <- Smooth_MSC_4dOb_2[rownames(Smooth_MSC_7dOb_2),]
Smooth_MSC_7dOb_3 <- Smooth_MSC_7dOb_3[rownames(Smooth_MSC_4dOb_3),]

Smooth_MSC_7dOb <- rbind(Smooth_MSC_7dOb_1, Smooth_MSC_7dOb_2, Smooth_MSC_7dOb_3) #5621
Smooth_MSC_4dOb <- rbind(Smooth_MSC_4dOb_1, Smooth_MSC_4dOb_2, Smooth_MSC_4dOb_3) #5621

Smooth_MSC_7dAd <- rbind(Smooth_MSC_7dAd_1, Smooth_MSC_7dAd_2, Smooth_MSC_7dAd_3) #3525
Smooth_MSC_4dAd <- rbind(Smooth_MSC_4dAd_1, Smooth_MSC_4dAd_2, Smooth_MSC_4dAd_3) #3525


# Split heatmap (artificially introduce white spaces)
h <- as.data.frame(matrix(0, nrow=200, ncol=200))
colnames(h) <- c(colnames(Smooth_MSC_7dAd), colnames(Smooth_MSC_4dAd))

#plot figure 3G - Ob
Heatmap(rbind(cbind(Smooth_MSC_7dOb_1, Smooth_MSC_4dOb_1), h, cbind(Smooth_MSC_7dOb_2, Smooth_MSC_4dOb_2), h, cbind(Smooth_MSC_7dOb_3, Smooth_MSC_4dOb_3)), cluster_columns=F, cluster_rows=F, column_title=paste("7dOb_vs_4dOb_no_genes_",nrow(Smooth_MSC_7dOb),"with_p_<",pvalue))

#plot figure 3J - Ad
Heatmap(rbind(cbind(Smooth_MSC_7dAd_1, Smooth_MSC_4dAd_1), h, cbind(Smooth_MSC_7dAd_2, Smooth_MSC_4dAd_2), h, cbind(Smooth_MSC_7dAd_3, Smooth_MSC_4dAd_3)), cluster_columns=F, cluster_rows=F, column_title=paste("7dAd_vs_4dAd_no_genes_",nrow(Smooth_MSC_7dAd),"with_p_<",pvalue))

### Figure 3H & 3K

# perform correlation across the heatmaps
# Ob - Figure 3H
row_corrs_1 <- mapply(function(x, y) cor(x, y), 
                    split(Smooth_MSC_7dOb_1, row(Smooth_MSC_7dOb_1)), 
                    split(Smooth_MSC_4dOb_1, row(Smooth_MSC_4dOb_1)))
row_corrs_2 <- mapply(function(x, y) cor(x, y), 
                      split(Smooth_MSC_7dOb_2, row(Smooth_MSC_7dOb_2)), 
                      split(Smooth_MSC_4dOb_2, row(Smooth_MSC_4dOb_2)))
row_corrs_3 <- mapply(function(x, y) cor(x, y), 
                      split(Smooth_MSC_7dOb_3, row(Smooth_MSC_7dOb_3)), 
                      split(Smooth_MSC_4dOb_3, row(Smooth_MSC_4dOb_3)))

row_corrs_1 <- row_corrs_1[order(row_corrs_1)]
row_corrs_2 <- row_corrs_2[order(row_corrs_2)]
row_corrs_3 <- row_corrs_3[order(row_corrs_3)]

plot(0,0,pch="",xlim=c(0,1),ylim=c(-1,1), xlab="",ylab="Pearson's correlation")
lines((1:length(row_corrs_1))/length(row_corrs_1),row_corrs_1, col="grey")
lines((1:length(row_corrs_2))/length(row_corrs_2),row_corrs_2, col="black")
lines((1:length(row_corrs_3))/length(row_corrs_3),row_corrs_3, col="green")

# Ad - Figure 3K
row_corrs_1 <- mapply(function(x, y) cor(x, y), 
                      split(Smooth_MSC_7dAd_1, row(Smooth_MSC_7dAd_1)), 
                      split(Smooth_MSC_4dAd_1, row(Smooth_MSC_4dAd_1)))
row_corrs_2 <- mapply(function(x, y) cor(x, y), 
                      split(Smooth_MSC_7dAd_2, row(Smooth_MSC_7dAd_2)), 
                      split(Smooth_MSC_4dAd_2, row(Smooth_MSC_4dAd_2)))
row_corrs_3 <- mapply(function(x, y) cor(x, y), 
                      split(Smooth_MSC_7dAd_3, row(Smooth_MSC_7dAd_3)), 
                      split(Smooth_MSC_4dAd_3, row(Smooth_MSC_4dAd_3)))

row_corrs_1 <- row_corrs_1[order(row_corrs_1)]
row_corrs_2 <- row_corrs_2[order(row_corrs_2)]
row_corrs_3 <- row_corrs_3[order(row_corrs_3)]

plot(0,0,pch="",xlim=c(0,1),ylim=c(-1,1), xlab="",ylab="Pearson's correlation")
lines((1:length(row_corrs_1))/length(row_corrs_1),row_corrs_1, col="grey")
lines((1:length(row_corrs_2))/length(row_corrs_2),row_corrs_2, col="black")
lines((1:length(row_corrs_3))/length(row_corrs_3),row_corrs_3, col="green")

# save Ob and Ad heatmaps for figure 3I and 3L
write.table(Smooth_MSC_7dOb_1, file="Smooth_MSC_7dOb_1.txt", col.names = T, row.names=T, quote=F, sep="\t")
write.table(Smooth_MSC_7dOb_2, file="Smooth_MSC_7dOb_2.txt", col.names = T, row.names=T, quote=F, sep="\t")
write.table(Smooth_MSC_7dOb_3, file="Smooth_MSC_7dOb_3.txt", col.names = T, row.names=T, quote=F, sep="\t")
write.table(Smooth_MSC_4dOb_1, file="Smooth_MSC_4dOb_1.txt", col.names = T, row.names=T, quote=F, sep="\t")
write.table(Smooth_MSC_4dOb_2, file="Smooth_MSC_4dOb_2.txt", col.names = T, row.names=T, quote=F, sep="\t")
write.table(Smooth_MSC_4dOb_3, file="Smooth_MSC_4dOb_3.txt", col.names = T, row.names=T, quote=F, sep="\t")
write.table(Smooth_MSC_7dAd_1, file="Smooth_MSC_7dAd_1.txt", col.names = T, row.names=T, quote=F, sep="\t")
write.table(Smooth_MSC_7dAd_2, file="Smooth_MSC_7dAd_2.txt", col.names = T, row.names=T, quote=F, sep="\t")
write.table(Smooth_MSC_7dAd_3, file="Smooth_MSC_7dAd_3.txt", col.names = T, row.names=T, quote=F, sep="\t")
write.table(Smooth_MSC_4dAd_1, file="Smooth_MSC_4dAd_1.txt", col.names = T, row.names=T, quote=F, sep="\t")
write.table(Smooth_MSC_4dAd_2, file="Smooth_MSC_4dAd_2.txt", col.names = T, row.names=T, quote=F, sep="\t")
write.table(Smooth_MSC_4dAd_3, file="Smooth_MSC_4dAd_3.txt", col.names = T, row.names=T, quote=F, sep="\t")

### Figure 3I, comparing scRNA-seq dynamics and bulk RNA-seq in PCA for Ob

# First bin the scRNA-seq dynamics from Figure 3G
Smooth_MSC_7dOb_1 <- read.delim("Smooth_MSC_7dOb_1.txt", h = T)
Smooth_MSC_7dOb_2 <- read.delim("Smooth_MSC_7dOb_2.txt", h = T)
Smooth_MSC_7dOb_3 <- read.delim("Smooth_MSC_7dOb_3.txt", h = T)

Smooth_MSC_7dOb <- rbind(Smooth_MSC_7dOb_1,Smooth_MSC_7dOb_2,Smooth_MSC_7dOb_3)

# use 10 bins
bin_size <- 20
n_bins <- ncol(Smooth_MSC_7dOb) / bin_size
# Use sapply to compute row means for each bin
binned_Smooth_MSC_7dOb <- sapply(1:n_bins, function(i) {
  cols <- ((i - 1) * bin_size + 1):(i * bin_size)
  rowMeans(Smooth_MSC_7dOb[, cols, drop = FALSE], na.rm = TRUE)
})
# add colnames to indicate matrix origin and bins
colnames(binned_Smooth_MSC_7dOb) <- paste("Smooth_7dOb_",1:5,sep="")

# do the same for 7dAd + 4d Ob
Smooth_MSC_4dOb_1 <- read.delim("Smooth_MSC_4dOb_1.txt", h = T)
Smooth_MSC_4dOb_2 <- read.delim("Smooth_MSC_4dOb_2.txt", h = T)
Smooth_MSC_4dOb_3 <- read.delim("Smooth_MSC_4dOb_3.txt", h = T)
Smooth_MSC_4dOb <- rbind(Smooth_MSC_4dOb_1,Smooth_MSC_4dOb_2,Smooth_MSC_4dOb_3)

# Use sapply to compute row means for each bin
binned_Smooth_MSC_4dOb <- sapply(1:n_bins, function(i) {
  cols <- ((i - 1) * bin_size + 1):(i * bin_size)
  rowMeans(Smooth_MSC_4dOb[, cols, drop = FALSE], na.rm = TRUE)
})

# add colnames to indicate matrix origin and bins
colnames(binned_Smooth_MSC_4dOb) <- paste("Smooth_4dOb_",1:5,sep="")

# Combine both binned scRNA-seq daynamics (7dOb and 7dAd + 4dOb)
Ob_Smooth <- cbind(binned_Smooth_MSC_7dOb,binned_Smooth_MSC_4dOb)

# Load bulk RNA-seq data
dds_mrg <- readRDS("dds_mrg.Rds")
RNA_tbl <- readRDS("RNA_tbl.Rds")

tmp <- data.frame(assay(rlogTransformation(dds_mrg, blind = FALSE)))
tmp$RefSeqID <- rownames(tmp)
tmp <- merge(tmp,RNA_tbl[,c("RefSeqID","Symbol")], by="RefSeqID")
rownames(tmp) <- tmp$Symbol
tmp$RefSeqID <- NULL
tmp$Symbol <- NULL

# select bulk RNA-seq samples of Ob lineage along with 7d Ad and switching 7dAd + 4dOb 
Ob <- tmp[,c(2:6,13:17,10,21,39,45)]

# create a helping vector that includes gene names which are in both data sets.
help <- rownames(Ob_Smooth)[rownames(Ob_Smooth) %in% rownames(tmp)]
Ob <- Ob[help,]
Ob_Smooth <- Ob_Smooth[help,]

# Generate PCA environment with bulk RNA-seq data 
pca <- prcomp(t(t(scale(t(Ob)))), scale.=F)

# center and scale bulk RNA-seq samples based on existing PCA
object2 <- scale(t(Ob_Smooth), center = pca$center, scale = pca$scale)
projected_scores2 <- object2 %*% pca$rotation

# Combine PCA scores
scores <- pca$x
scores <- rbind(scores,projected_scores2)
# Collect PCA scores and variance for plotting 
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
dat <- data.frame(PC1=scores[,1], PC2=scores[,2], PC3=scores[,3], PC4=scores[,4], name=c(colnames(Ob),colnames(Ob_Smooth)))
attr(dat, "percentVar") <- percentVar[1:4]

# Make extra column for sample and time point
dat$Group <- c('Ob','Ob','Ob','Ob','Msc','Ob','Ob','Ob','Ob','Msc','Ad','Ad','Ob','Ob','Smooth_7dOb','Smooth_7dOb','Smooth_7dOb','Smooth_7dOb','Smooth_7dOb','Smooth_4dOb','Smooth_4dOb','Smooth_4dOb','Smooth_4dOb','Smooth_4dOb')
dat$Timepoints <- as.factor(c(5,4,3,2,1,5,4,3,2,1,6,6,7,7,1,2,3,4,5,1,2,3,4,5))

# scale axis to get a squared graph
ratio.values <- (max(dat$PC1)-min(dat$PC1))/(max(dat$PC2)-min(dat$PC2))
#plot results
library(ggplot2)
ggplot(data=dat, aes_string(x="PC1", y="PC2", shape= "Group",color="Timepoints"))+ geom_point(size=3) + coord_fixed(ratio.values) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance"))



### Figure 3L, comparing scRNA-seq dynamics and bulk RNA-seq in PCA for Ad

# First bin the scRNA-seq dynamics from Figure 3G
Smooth_MSC_7dAd_1 <- read.delim("Smooth_MSC_7dAd_1.txt", h = T)
Smooth_MSC_7dAd_2 <- read.delim("Smooth_MSC_7dAd_2.txt", h = T)
Smooth_MSC_7dAd_3 <- read.delim("Smooth_MSC_7dAd_3.txt", h = T)

Smooth_MSC_7dAd <- rbind(Smooth_MSC_7dAd_1,Smooth_MSC_7dAd_2,Smooth_MSC_7dAd_3)

# use 10 bins
bin_size <- 20
n_bins <- ncol(Smooth_MSC_7dAd) / bin_size
# Use sapply to compute row means for each bin
binned_Smooth_MSC_7dAd <- sapply(1:n_bins, function(i) {
  cols <- ((i - 1) * bin_size + 1):(i * bin_size)
  rowMeans(Smooth_MSC_7dAd[, cols, drop = FALSE], na.rm = TRUE)
})
# add colnames to indicate matrix origin and bins
colnames(binned_Smooth_MSC_7dAd) <- paste("Smooth_7dOb_",1:5,sep="")

# do the same for 7dAd + 4d Ob
Smooth_MSC_4dAd_1 <- read.delim("Smooth_MSC_4dAd_1.txt", h = T)
Smooth_MSC_4dAd_2 <- read.delim("Smooth_MSC_4dAd_2.txt", h = T)
Smooth_MSC_4dAd_3 <- read.delim("Smooth_MSC_4dAd_3.txt", h = T)
Smooth_MSC_4dAd <- rbind(Smooth_MSC_4dAd_1,Smooth_MSC_4dAd_2,Smooth_MSC_4dAd_3)

# Use sapply to compute row means for each bin
binned_Smooth_MSC_4dAd <- sapply(1:n_bins, function(i) {
  cols <- ((i - 1) * bin_size + 1):(i * bin_size)
  rowMeans(Smooth_MSC_4dAd[, cols, drop = FALSE], na.rm = TRUE)
})

# add colnames to indicate matrix origin and bins
colnames(binned_Smooth_MSC_4dAd) <- paste("Smooth_4dOb_",1:5,sep="")

# Combine both binned scRNA-seq daynamics (7dOb and 7dAd + 4dOb)
Ad_Smooth <- cbind(binned_Smooth_MSC_7dAd,binned_Smooth_MSC_4dAd)

# Load bulk RNA-seq data
dds_mrg <- readRDS("dds_mrg.Rds")
RNA_tbl <- readRDS("RNA_tbl.Rds")

tmp <- data.frame(assay(rlogTransformation(dds_mrg, blind = FALSE)))
tmp$RefSeqID <- rownames(tmp)
tmp <- merge(tmp,RNA_tbl[,c("RefSeqID","Symbol")], by="RefSeqID")
rownames(tmp) <- tmp$Symbol
tmp$RefSeqID <- NULL
tmp$Symbol <- NULL

# select bulk RNA-seq samples of Ob lineage along with 7d Ad and switching 7dAd + 4dOb 
Ad <- tmp[,c(6:10,17:21,2,13,38,44)]

# create a helping vector that includes gene names which are in both data sets.
help <- rownames(Ad_Smooth)[rownames(Ad_Smooth) %in% rownames(tmp)]
Ad <- Ad[help,]
Ad_Smooth <- Ad_Smooth[help,]

# Generate PCA environment with bulk RNA-seq data 
pca <- prcomp(t(t(scale(t(Ad)))), scale.=F)

# center and scale scRNA-seq samples based on existing PCA
object2 <- scale(t(Ad_Smooth), center = pca$center, scale = pca$scale)
projected_scores2 <- object2 %*% pca$rotation

# Combine PCA scores
scores <- pca$x
scores <- rbind(scores,projected_scores2)
# Collect PCA scores and variance for plotting 
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
dat <- data.frame(PC1=scores[,1], PC2=scores[,2], PC3=scores[,3], PC4=scores[,4], name=c(colnames(Ad),colnames(Ad_Smooth)))
attr(dat, "percentVar") <- percentVar[1:4]

# Make extra column for sample and time point
dat$Group <- c('Msc','Ad','Ad','Ad','Ad','Msc','Ad','Ad','Ad','Ad','Ob','Ob','Ad','Ad','Smooth_7dAd','Smooth_7dAd','Smooth_7dAd','Smooth_7dAd','Smooth_7dAd','Smooth_4dAd','Smooth_4dAd','Smooth_4dAd','Smooth_4dAd','Smooth_4dAd')
dat$Timepoints <- as.factor(c(1,2,3,4,5,1,2,3,4,5,6,6,7,7,1,2,3,4,5,1,2,3,4,5))

# scale axis to get a squared graph
ratio.values <- (max(dat$PC1)-min(dat$PC1))/(max(dat$PC2)-min(dat$PC2))
#plot results
library(ggplot2)
ggplot(data=dat, aes_string(x="PC1", y="PC2", shape= "Group",color="Timepoints"))+ geom_point(size=3) + coord_fixed(ratio.values) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance"))

