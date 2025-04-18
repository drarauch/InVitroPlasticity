### Figure 3
# load the following objects (see DataPreparation.R for more information on how these were generated)
seurat_all <- readRDS("sc_seurat_diff_scales.Rds")

### Figure 3B
# UMAP-plot of seurat object
DimPlot(seurat_all, group.by="Dataset")

### Figure 3C
# Violin plot of Osteoblast (ALPL & LMCD1) and Adipocyte (PLIN1 & ADIPOQ) markers within each dataset
VlnPlot(seurat_all, c("ALPL", "LMCD1"), group.by="Dataset")
VlnPlot(seurat_all, c("PLIN1","ADIPOQ"), group.by="Dataset")

### Figure 3D
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
#saveRDS(fitGAM_7dAd, "x7dAd_only_woProf_fitGAM_2")
#saveRDS(fitGAM_4dAd, "x7dOb4dAd_only_woProf_fitGAM_2")
#saveRDS(fitGAM_7dOb, "x7dOb_sce_fitGAM_new")
#saveRDS(fitGAM_4dOb, "x7dAd4dOb_only_woProf_fitGAM_new")

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

## load objects
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

#Subset Ass_xx to regulated common genes
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

# Average across replicates and scale
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

# Seriate the results
Smooth_MSC_7dAd_1 <- Smooth_MSC_7dAd_1[ get_order(seriate(Smooth_MSC_7dAd_1, method="PCA_angle")),]
Smooth_MSC_7dAd_2 <- Smooth_MSC_7dAd_2[ get_order(seriate(Smooth_MSC_7dAd_2, method="PCA_angle")),]
Smooth_MSC_7dAd_3 <- Smooth_MSC_7dAd_3[ get_order(seriate(Smooth_MSC_7dAd_3, method="PCA_angle")),]

Smooth_MSC_4dAd_1 <- Smooth_MSC_4dAd_1[ get_order(seriate(Smooth_MSC_4dAd_1, method="PCA_angle")),]
Smooth_MSC_4dAd_2 <- Smooth_MSC_4dAd_2[ get_order(seriate(Smooth_MSC_4dAd_2, method="PCA_angle")),]
Smooth_MSC_4dAd_3 <- Smooth_MSC_4dAd_3[ get_order(seriate(Smooth_MSC_4dAd_3, method="PCA_angle")),]

Smooth_MSC_7dOb_1 <- Smooth_MSC_7dOb_1[ get_order(seriate(Smooth_MSC_7dOb_1, method="PCA_angle")),]
Smooth_MSC_7dOb_2 <- Smooth_MSC_7dOb_2[ get_order(seriate(Smooth_MSC_7dOb_2, method="PCA_angle")),]
Smooth_MSC_7dOb_3 <- Smooth_MSC_7dOb_3[ get_order(seriate(Smooth_MSC_7dOb_3, method="PCA_angle")),]

Smooth_MSC_4dOb_1 <- Smooth_MSC_4dOb_1[ get_order(seriate(Smooth_MSC_4dOb_1, method="PCA_angle")),]
Smooth_MSC_4dOb_2 <- Smooth_MSC_4dOb_2[ get_order(seriate(Smooth_MSC_4dOb_2, method="PCA_angle")),]
Smooth_MSC_4dOb_3 <- Smooth_MSC_4dOb_3[ get_order(seriate(Smooth_MSC_4dOb_3, method="PCA_angle")),]

Smooth_MSC_7dAd_1 <- as.data.frame(Smooth_MSC_7dAd_1)
Smooth_MSC_7dAd_2 <- as.data.frame(Smooth_MSC_7dAd_2)
Smooth_MSC_7dAd_3 <- as.data.frame(Smooth_MSC_7dAd_3)

Smooth_MSC_4dAd_1 <- as.data.frame(Smooth_MSC_4dAd_1)
Smooth_MSC_4dAd_2 <- as.data.frame(Smooth_MSC_4dAd_2)
Smooth_MSC_4dAd_3 <- as.data.frame(Smooth_MSC_4dAd_3)

Smooth_MSC_7dOb_1 <- as.data.frame(Smooth_MSC_7dOb_1)
Smooth_MSC_7dOb_2 <- as.data.frame(Smooth_MSC_7dOb_2)
Smooth_MSC_7dOb_3 <- as.data.frame(Smooth_MSC_7dOb_3)

Smooth_MSC_4dOb_1 <- as.data.frame(Smooth_MSC_4dOb_1)
Smooth_MSC_4dOb_2 <- as.data.frame(Smooth_MSC_4dOb_2)
Smooth_MSC_4dOb_3 <- as.data.frame(Smooth_MSC_4dOb_3)

Smooth_MSC_7dAd_1$genes <- rownames(Smooth_MSC_7dAd_1)
Smooth_MSC_7dAd_2$genes <- rownames(Smooth_MSC_7dAd_2)
Smooth_MSC_7dAd_3$genes <- rownames(Smooth_MSC_7dAd_3)

Smooth_MSC_4dAd_1$genes <- rownames(Smooth_MSC_4dAd_1)
Smooth_MSC_4dAd_2$genes <- rownames(Smooth_MSC_4dAd_2)
Smooth_MSC_4dAd_3$genes <- rownames(Smooth_MSC_4dAd_3)

Smooth_MSC_7dOb_1$genes <- rownames(Smooth_MSC_7dOb_1)
Smooth_MSC_7dOb_2$genes <- rownames(Smooth_MSC_7dOb_2)
Smooth_MSC_7dOb_3$genes <- rownames(Smooth_MSC_7dOb_3)

Smooth_MSC_4dOb_1$genes <- rownames(Smooth_MSC_4dOb_1)
Smooth_MSC_4dOb_2$genes <- rownames(Smooth_MSC_4dOb_2)
Smooth_MSC_4dOb_3$genes <- rownames(Smooth_MSC_4dOb_3)

Smooth_MSC_4dAd_1 <- Smooth_MSC_4dAd_1  %>% slice(match(as.matrix(Smooth_MSC_7dAd_1), genes))
Smooth_MSC_4dAd_2 <- Smooth_MSC_4dAd_2  %>% slice(match(as.matrix(Smooth_MSC_7dAd_2), genes))
Smooth_MSC_7dAd_3 <- Smooth_MSC_7dAd_3  %>% slice(match(as.matrix(Smooth_MSC_4dAd_3), genes))

Smooth_MSC_4dOb_1 <- Smooth_MSC_4dOb_1  %>% slice(match(as.matrix(Smooth_MSC_7dOb_1), genes))
Smooth_MSC_4dOb_2 <- Smooth_MSC_4dOb_2  %>% slice(match(as.matrix(Smooth_MSC_7dOb_2), genes))
Smooth_MSC_7dOb_3 <- Smooth_MSC_7dOb_3  %>% slice(match(as.matrix(Smooth_MSC_4dOb_3), genes))

Smooth_MSC_4dAd_1$genes <- NULL
Smooth_MSC_4dAd_2$genes <- NULL
Smooth_MSC_4dAd_3$genes <- NULL

Smooth_MSC_7dAd_1$genes <- NULL
Smooth_MSC_7dAd_2$genes <- NULL
Smooth_MSC_7dAd_3$genes <- NULL

Smooth_MSC_4dOb_1$genes <- NULL
Smooth_MSC_4dOb_2$genes <- NULL
Smooth_MSC_4dOb_3$genes <- NULL

Smooth_MSC_7dOb_1$genes <- NULL
Smooth_MSC_7dOb_2$genes <- NULL
Smooth_MSC_7dOb_3$genes <- NULL

Smooth_MSC_7dOb <- rbind(Smooth_MSC_7dOb_1, Smooth_MSC_7dOb_2, Smooth_MSC_7dOb_3) #5621
Smooth_MSC_4dOb <- rbind(Smooth_MSC_4dOb_1, Smooth_MSC_4dOb_2, Smooth_MSC_4dOb_3) #5621

Smooth_MSC_7dAd <- rbind(Smooth_MSC_7dAd_1, Smooth_MSC_7dAd_2, Smooth_MSC_7dAd_3) #3525
Smooth_MSC_4dAd <- rbind(Smooth_MSC_4dAd_1, Smooth_MSC_4dAd_2, Smooth_MSC_4dAd_3) #3525


# Split heatmap
h <- as.data.frame(matrix(0, nrow=200, ncol=200))
colnames(h) <- c(colnames(Smooth_MSC_7dAd), colnames(Smooth_MSC_4dAd))

#plot figure 3G
Heatmap(rbind(cbind(Smooth_MSC_7dOb_1, Smooth_MSC_4dOb_1), h, cbind(Smooth_MSC_7dOb_2, Smooth_MSC_4dOb_2), h, cbind(Smooth_MSC_7dOb_3, Smooth_MSC_4dOb_3)), cluster_columns=F, cluster_rows=F, column_title=paste("7dOb_vs_4dOb_no_genes_",nrow(Smooth_MSC_7dOb),"with_p_<",pvalue))

#plot figure 3J
Heatmap(rbind(cbind(Smooth_MSC_7dAd_1, Smooth_MSC_4dAd_1), h, cbind(Smooth_MSC_7dAd_2, Smooth_MSC_4dAd_2), h, cbind(Smooth_MSC_7dAd_3, Smooth_MSC_4dAd_3)), cluster_columns=F, cluster_rows=F, column_title=paste("7dAd_vs_4dAd_no_genes_",nrow(Smooth_MSC_7dAd),"with_p_<",pvalue))

### Figure 3H & 3K
#order genes by pval
Ass_7dAd_1 <- Ass_7dAd_1[order(Ass_7dAd_1$pvalue),]
Ass_7dAd_2 <- Ass_7dAd_2[order(Ass_7dAd_2$pvalue),]
Ass_7dAd_3 <- Ass_7dAd_3[order(Ass_7dAd_3$pvalue),]
Ass_4dAd_1 <- Ass_4dAd_1[order(Ass_4dAd_1$pvalue),]
Ass_4dAd_2 <- Ass_4dAd_2[order(Ass_4dAd_2$pvalue),]
Ass_4dAd_3 <- Ass_4dAd_3[order(Ass_4dAd_3$pvalue),]
Ass_7dOb_1 <- Ass_7dOb_1[order(Ass_7dOb_1$pvalue),]
Ass_7dOb_2 <- Ass_7dOb_2[order(Ass_7dOb_2$pvalue),]
Ass_7dOb_3 <- Ass_7dOb_3[order(Ass_7dOb_3$pvalue),]
Ass_4dOb_1 <- Ass_4dOb_1[order(Ass_4dOb_1$pvalue),]
Ass_4dOb_2 <- Ass_4dOb_2[order(Ass_4dOb_2$pvalue),]
Ass_4dOb_3 <- Ass_4dOb_3[order(Ass_4dOb_3$pvalue),]

#add rank values
Ass_7dAd_1$rank_pval <- c(1:nrow(Ass_7dAd_1))
Ass_7dAd_2$rank_pval <- c(1:nrow(Ass_7dAd_2))
Ass_7dAd_3$rank_pval <- c(1:nrow(Ass_7dAd_3))
Ass_4dAd_1$rank_pval <- c(1:nrow(Ass_4dAd_1))
Ass_4dAd_2$rank_pval <- c(1:nrow(Ass_4dAd_2))
Ass_4dAd_3$rank_pval <- c(1:nrow(Ass_4dAd_3))
Ass_7dOb_1$rank_pval <- c(1:nrow(Ass_4dOb_1))
Ass_7dOb_2$rank_pval <- c(1:nrow(Ass_4dOb_2))
Ass_7dOb_3$rank_pval <- c(1:nrow(Ass_4dOb_3))
Ass_4dOb_1$rank_pval <- c(1:nrow(Ass_4dOb_1))
Ass_4dOb_2$rank_pval <- c(1:nrow(Ass_4dOb_2))
Ass_4dOb_3$rank_pval <- c(1:nrow(Ass_4dOb_3))

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

# Average across replicates and scale
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

# Seriate the results
Smooth_MSC_7dAd_1 <- Smooth_MSC_7dAd_1[ get_order(seriate(Smooth_MSC_7dAd_1, method="PCA_angle")),]
Smooth_MSC_7dAd_2 <- Smooth_MSC_7dAd_2[ get_order(seriate(Smooth_MSC_7dAd_2, method="PCA_angle")),]
Smooth_MSC_7dAd_3 <- Smooth_MSC_7dAd_3[ get_order(seriate(Smooth_MSC_7dAd_3, method="PCA_angle")),]
Smooth_MSC_4dAd_1 <- Smooth_MSC_4dAd_1[ get_order(seriate(Smooth_MSC_4dAd_1, method="PCA_angle")),]
Smooth_MSC_4dAd_2 <- Smooth_MSC_4dAd_2[ get_order(seriate(Smooth_MSC_4dAd_2, method="PCA_angle")),]
Smooth_MSC_4dAd_3 <- Smooth_MSC_4dAd_3[ get_order(seriate(Smooth_MSC_4dAd_3, method="PCA_angle")),]
Smooth_MSC_7dOb_1 <- Smooth_MSC_7dOb_1[ get_order(seriate(Smooth_MSC_7dOb_1, method="PCA_angle")),]
Smooth_MSC_7dOb_2 <- Smooth_MSC_7dOb_2[ get_order(seriate(Smooth_MSC_7dOb_2, method="PCA_angle")),]
Smooth_MSC_7dOb_3 <- Smooth_MSC_7dOb_3[ get_order(seriate(Smooth_MSC_7dOb_3, method="PCA_angle")),]
Smooth_MSC_4dOb_1 <- Smooth_MSC_4dOb_1[ get_order(seriate(Smooth_MSC_4dOb_1, method="PCA_angle")),]
Smooth_MSC_4dOb_2 <- Smooth_MSC_4dOb_2[ get_order(seriate(Smooth_MSC_4dOb_2, method="PCA_angle")),]
Smooth_MSC_4dOb_3 <- Smooth_MSC_4dOb_3[ get_order(seriate(Smooth_MSC_4dOb_3, method="PCA_angle")),]

Smooth_MSC_7dAd_1 <- as.data.frame(Smooth_MSC_7dAd_1)
Smooth_MSC_7dAd_2 <- as.data.frame(Smooth_MSC_7dAd_2)
Smooth_MSC_7dAd_3 <- as.data.frame(Smooth_MSC_7dAd_3)
Smooth_MSC_4dAd_1 <- as.data.frame(Smooth_MSC_4dAd_1)
Smooth_MSC_4dAd_2 <- as.data.frame(Smooth_MSC_4dAd_2)
Smooth_MSC_4dAd_3 <- as.data.frame(Smooth_MSC_4dAd_3)
Smooth_MSC_7dOb_1 <- as.data.frame(Smooth_MSC_7dOb_1)
Smooth_MSC_7dOb_2 <- as.data.frame(Smooth_MSC_7dOb_2)
Smooth_MSC_7dOb_3 <- as.data.frame(Smooth_MSC_7dOb_3)
Smooth_MSC_4dOb_1 <- as.data.frame(Smooth_MSC_4dOb_1)
Smooth_MSC_4dOb_2 <- as.data.frame(Smooth_MSC_4dOb_2)
Smooth_MSC_4dOb_3 <- as.data.frame(Smooth_MSC_4dOb_3)

Smooth_MSC_7dAd_1$genes <- rownames(Smooth_MSC_7dAd_1)
Smooth_MSC_7dAd_2$genes <- rownames(Smooth_MSC_7dAd_2)
Smooth_MSC_7dAd_3$genes <- rownames(Smooth_MSC_7dAd_3)
Smooth_MSC_4dAd_1$genes <- rownames(Smooth_MSC_4dAd_1)
Smooth_MSC_4dAd_2$genes <- rownames(Smooth_MSC_4dAd_2)
Smooth_MSC_4dAd_3$genes <- rownames(Smooth_MSC_4dAd_3)
Smooth_MSC_7dOb_1$genes <- rownames(Smooth_MSC_7dOb_1)
Smooth_MSC_7dOb_2$genes <- rownames(Smooth_MSC_7dOb_2)
Smooth_MSC_7dOb_3$genes <- rownames(Smooth_MSC_7dOb_3)
Smooth_MSC_4dOb_1$genes <- rownames(Smooth_MSC_4dOb_1)
Smooth_MSC_4dOb_2$genes <- rownames(Smooth_MSC_4dOb_2)
Smooth_MSC_4dOb_3$genes <- rownames(Smooth_MSC_4dOb_3)

Smooth_MSC_4dAd_1 <- Smooth_MSC_4dAd_1  %>% slice(match(as.matrix(Smooth_MSC_7dAd_1), genes))
Smooth_MSC_4dAd_2 <- Smooth_MSC_4dAd_2  %>% slice(match(as.matrix(Smooth_MSC_7dAd_2), genes))
Smooth_MSC_7dAd_3 <- Smooth_MSC_7dAd_3  %>% slice(match(as.matrix(Smooth_MSC_4dAd_3), genes))
Smooth_MSC_4dOb_1 <- Smooth_MSC_4dOb_1  %>% slice(match(as.matrix(Smooth_MSC_7dOb_1), genes))
Smooth_MSC_4dOb_2 <- Smooth_MSC_4dOb_2  %>% slice(match(as.matrix(Smooth_MSC_7dOb_2), genes))
Smooth_MSC_7dOb_3 <- Smooth_MSC_7dOb_3  %>% slice(match(as.matrix(Smooth_MSC_4dOb_3), genes))

Smooth_MSC_4dAd_1$genes <- NULL
Smooth_MSC_4dAd_2$genes <- NULL
Smooth_MSC_4dAd_3$genes <- NULL
Smooth_MSC_7dAd_1$genes <- NULL
Smooth_MSC_7dAd_2$genes <- NULL
Smooth_MSC_7dAd_3$genes <- NULL
Smooth_MSC_4dOb_1$genes <- NULL
Smooth_MSC_4dOb_2$genes <- NULL
Smooth_MSC_4dOb_3$genes <- NULL
Smooth_MSC_7dOb_1$genes <- NULL
Smooth_MSC_7dOb_2$genes <- NULL
Smooth_MSC_7dOb_3$genes <- NULL

Smooth_MSC_4dAd_1 <- t(Smooth_MSC_4dAd_1)
Smooth_MSC_4dAd_2 <- t(Smooth_MSC_4dAd_2)
Smooth_MSC_4dAd_3 <- t(Smooth_MSC_4dAd_3)
Smooth_MSC_7dAd_1 <- t(Smooth_MSC_7dAd_1)
Smooth_MSC_7dAd_2 <- t(Smooth_MSC_7dAd_2)
Smooth_MSC_7dAd_3 <- t(Smooth_MSC_7dAd_3)
Smooth_MSC_4dOb_1 <- t(Smooth_MSC_4dOb_1)
Smooth_MSC_4dOb_2 <- t(Smooth_MSC_4dOb_2)
Smooth_MSC_4dOb_3 <- t(Smooth_MSC_4dOb_3)
Smooth_MSC_7dOb_1 <- t(Smooth_MSC_7dOb_1)
Smooth_MSC_7dOb_2 <- t(Smooth_MSC_7dOb_2)
Smooth_MSC_7dOb_3 <- t(Smooth_MSC_7dOb_3)

#7dOb & 7dAd 4dOb
cor_list_unfiltered_Ob_1 <- cor(Smooth_MSC_7dOb_1, Smooth_MSC_4dOb_1)
cor_list_unfiltered_Ob_2 <- cor(Smooth_MSC_7dOb_2, Smooth_MSC_4dOb_2)
cor_list_unfiltered_Ob_3 <- cor(Smooth_MSC_7dOb_3, Smooth_MSC_4dOb_3)

cor_list_Ob_1 <- NA
cor_list_Ob_2 <- NA
cor_list_Ob_3 <- NA

cor_list_Ob_1 <- as.data.frame(cor_list_Ob_1)
cor_list_Ob_2 <- as.data.frame(cor_list_Ob_2)
cor_list_Ob_3 <- as.data.frame(cor_list_Ob_3)

for (i in 1:nrow(cor_list_unfiltered_Ob_1)) {
  
  cor_list_Ob_1[i,] <- cor_list_unfiltered_Ob_1[i,i]
}
for (i in 1:nrow(cor_list_unfiltered_Ob_2)) {
  
  cor_list_Ob_2[i,] <- cor_list_unfiltered_Ob_2[i,i]
}
for (i in 1:nrow(cor_list_unfiltered_Ob_3)) {
  
  cor_list_Ob_3[i,] <- cor_list_unfiltered_Ob_3[i,i]
}

rownames(cor_list_Ob_1) <- colnames(Smooth_MSC_7dOb_1)
rownames(cor_list_Ob_2) <- colnames(Smooth_MSC_7dOb_2)
rownames(cor_list_Ob_3) <- colnames(Smooth_MSC_7dOb_3)

colnames(cor_list_Ob_1) <- "cor"
colnames(cor_list_Ob_2) <- "cor"
colnames(cor_list_Ob_3) <- "cor"
cor_list_Ob_1$genes <- rownames(cor_list_Ob_1)
cor_list_Ob_2$genes <- rownames(cor_list_Ob_2)
cor_list_Ob_3$genes <- rownames(cor_list_Ob_3)
cor_list_Ob_1 <- as.data.frame(cor_list_Ob_1[order(cor_list_Ob_1$cor, decreasing=TRUE),])
cor_list_Ob_2 <- as.data.frame(cor_list_Ob_2[order(cor_list_Ob_2$cor, decreasing=TRUE),])
cor_list_Ob_3 <- as.data.frame(cor_list_Ob_3[order(cor_list_Ob_3$cor, decreasing=TRUE),])
cor_list_Ob_1$rank <- c(nrow(cor_list_Ob_1):1)
cor_list_Ob_2$rank <- c(nrow(cor_list_Ob_2):1)
cor_list_Ob_3$rank <- c(nrow(cor_list_Ob_3):1)

#Plot figure 3H
#Grey
plot(cor_list_Ob_1$rank/nrow(cor_list_Ob_1), cor_list_Ob_1$cor, main = c("no. genes",nrow(cor_list_Ob_1)), ylim = c(-1,1), xlim=c(0,1))
#Black
plot(cor_list_Ob_2$rank/nrow(cor_list_Ob_2), cor_list_Ob_2$cor, main = c("no. genes",nrow(cor_list_Ob_2)), ylim = c(-1,1), xlim=c(0,1))
#Green
plot(cor_list_Ob_3$rank/nrow(cor_list_Ob_3), cor_list_Ob_3$cor, main = c("no. genes",nrow(cor_list_Ob_3)), ylim = c(-1,1), xlim=c(0,1))

#7dAd & 7dOb 4dAd
cor_list_unfiltered_Ad_1 <- cor(Smooth_MSC_7dAd_1, Smooth_MSC_4dAd_1)
cor_list_unfiltered_Ad_2 <- cor(Smooth_MSC_7dAd_2, Smooth_MSC_4dAd_2)
cor_list_unfiltered_Ad_3 <- cor(Smooth_MSC_7dAd_3, Smooth_MSC_4dAd_3)

cor_list_Ad_1 <- NA
cor_list_Ad_2 <- NA
cor_list_Ad_3 <- NA

cor_list_Ad_1 <- as.data.frame(cor_list_Ad_1)
cor_list_Ad_2 <- as.data.frame(cor_list_Ad_2)
cor_list_Ad_3 <- as.data.frame(cor_list_Ad_3)

for (i in 1:nrow(cor_list_unfiltered_Ad_1)) {
  
  cor_list_Ad_1[i,] <- cor_list_unfiltered_Ad_1[i,i]
}
for (i in 1:nrow(cor_list_unfiltered_Ad_2)) {
  
  cor_list_Ad_2[i,] <- cor_list_unfiltered_Ad_2[i,i]
}
for (i in 1:nrow(cor_list_unfiltered_Ad_3)) {
  
  cor_list_Ad_3[i,] <- cor_list_unfiltered_Ad_3[i,i]
}

rownames(cor_list_Ad_1) <- colnames(Smooth_MSC_7dAd_1)
rownames(cor_list_Ad_2) <- colnames(Smooth_MSC_7dAd_2)
rownames(cor_list_Ad_3) <- colnames(Smooth_MSC_7dAd_3)

colnames(cor_list_Ad_1) <- "cor"
colnames(cor_list_Ad_2) <- "cor"
colnames(cor_list_Ad_3) <- "cor"
cor_list_Ad_1$genes <- rownames(cor_list_Ad_1)
cor_list_Ad_2$genes <- rownames(cor_list_Ad_2)
cor_list_Ad_3$genes <- rownames(cor_list_Ad_3)
cor_list_Ad_1 <- as.data.frame(cor_list_Ad_1[order(cor_list_Ad_1$cor, decreasing=TRUE),])
cor_list_Ad_2 <- as.data.frame(cor_list_Ad_2[order(cor_list_Ad_2$cor, decreasing=TRUE),])
cor_list_Ad_3 <- as.data.frame(cor_list_Ad_3[order(cor_list_Ad_3$cor, decreasing=TRUE),])
cor_list_Ad_1$rank <- c(nrow(cor_list_Ad_1):1)
cor_list_Ad_2$rank <- c(nrow(cor_list_Ad_2):1)
cor_list_Ad_3$rank <- c(nrow(cor_list_Ad_3):1)

#Plot figure 3K
#Grey
plot(cor_list_Ad_1$rank/nrow(cor_list_Ad_1), cor_list_Ad_1$cor, main = c("no. genes",nrow(cor_list_Ad_1)), ylim = c(-1,1), xlim=c(0,1))
#Black
plot(cor_list_Ad_2$rank/nrow(cor_list_Ad_2), cor_list_Ad_2$cor, main = c("no. genes",nrow(cor_list_Ad_2)), ylim = c(-1,1), xlim=c(0,1))
#Green
plot(cor_list_Ad_3$rank/nrow(cor_list_Ad_3), cor_list_Ad_3$cor, main = c("no. genes",nrow(cor_list_Ad_3)), ylim = c(-1,1), xlim=c(0,1))

### Figure 3I & 3L

dds_mrg <- readRDS("dds_mrg.Rds")
RNA_tbl <- readRDS("RNA_tbl3.Rds")

tmp <- data.frame(assay(rlogTransformation(dds_mrg, blind = FALSE)))

#7dOb + 7dAd4dOb
tmp <- tmp[, c("Msc_1", "Msc_2", "Msc_3","X4hOb_1", "X4hOb_2", "X4hOb_3", "X1dOb_1", "X1dOb_2", "X1dOb_3", "X3dOb_1", "X3dOb_2", "X3dOb_3", "X7dOb_1", "X7dOb_2", "X7dOb_3", "X14dOb_1", "X14dOb_2", "X14dOb_3", "X7dAd_1", "X7dAd_2", "X7dAd_3", "X7Ad_4dOb_1", "X7dAd_4dOb_2")]
tmp <- na.omit(tmp)

#Get average values from multiple replicates
a <- 1:3
tmp$Msc <- NA
for (i in 1:nrow(tmp)) 
  
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*0]))
}
tmp$X4hOb <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*1]))
}
tmp$X1dOb <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*2]))
}
tmp$X3dOb <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*3]))
}
tmp$X7dOb <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*4]))
}
tmp$X14dOb <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*5]))
}
tmp$X7dAd <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*6]))
}
tmp$X7dAd_4dOb <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,22:23]))
}

tmp_cond <- tmp[,24:ncol(tmp)] 
tmp_cond$genes <- rownames(tmp_cond)

#Subset to genes from single cell object
x <- readRDS("A2.Rds")
colnames(x)[2] <- "genes"
tmptmp <- tmp_cond
tmp2 <- merge(tmptmp, x, by="genes")
tmp2$genes <- NULL
colnames(tmp2)[9] <- "genes"

#plot figure 3I
#x <- readRDS("genes_7dOb_only") #use for upper part 
#x <- readRDS("genes_Ob_common") #use for middle part
x <- readRDS("genes_4dOb_only") #use for lower part

tmp3 <- merge(tmp2, x, by="genes")
tmp3 <- tmp3  %>% slice(match(x$genes, genes))
rownames(tmp3) <- tmp3$genes
tmp3$genes <- NULL

#scale data
tmp3 <- t(scale(t(tmp3))) 

#set limits
tmp3[ tmp3 > 2] <- 2
tmp3[ tmp3 < -2] <- -2

pheatmap(tmp3, cluster_cols = F, cluster_rows = F, col=rev(myCol))

#7dAd + 7dOb4dAd
tmp <- tmp[, c("Msc_1", "Msc_2", "Msc_3", "X4hAd_1", "X4hAd_2", "X4hAd_3", "X1dAd_1", "X1dAd_2", "X1dAd_3", "X3dAd_1", "X3dAd_2", "X3dAd_3", "X7dAd_1", "X7dAd_2", "X7dAd_3", "X14dAd_1", "X14dAd_2", "X14dAd_3", "X7dOb_1", "X7dOb_2", "X7dOb_3", "X7dOb_4dAd_1", "X7dOb_4dAd_2")]
tmp <- na.omit(tmp)

#Get average values from multiple replicates
a <- 1:3
tmp$Msc <- NA
for (i in 1:nrow(tmp)) 
  
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*0]))
}
tmp$X4hAd <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*1]))
}
tmp$X1dAd <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*2]))
}
tmp$X3dAd <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*3]))
}
tmp$X7dAd <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*4]))
}
tmp$X14dAd <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*5]))
}
tmp$X7dOb <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,a+3*6]))
}
tmp$X7dOb_4dAd <- NA
for (i in 1:nrow(tmp)) 
{ 
  tmp[i,ncol(tmp)] <- mean(as.numeric(tmp[i,22:23]))
}

tmp_cond <- tmp[,24:ncol(tmp)] 
tmp_cond$genes <- rownames(tmp_cond)

#Subset to genes from single cell object
x <- readRDS("A2.Rds")
colnames(x)[2] <- "genes"
tmptmp <- tmp_cond
tmp2 <- merge(tmptmp, x, by="genes")
tmp2$genes <- NULL
colnames(tmp2)[9] <- "genes"

#plot figure 3L
x <- readRDS("genes_7dAd_only") #use for upper part 
x <- readRDS("genes_Ad_common") #use for middle part
x <- readRDS("genes_4dAd_only") #use for lower part

tmp3 <- merge(tmp2, x, by="genes")
tmp3 <- tmp3  %>% slice(match(x$genes, genes))
rownames(tmp3) <- tmp3$genes
tmp3$genes <- NULL

#scale data
tmp3 <- t(scale(t(tmp3))) 

#set limits
tmp3[ tmp3 > 2] <- 2
tmp3[ tmp3 < -2] <- -2

pheatmap(tmp3, cluster_cols = F, cluster_rows = F, col=rev(myCol))

