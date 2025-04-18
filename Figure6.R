### Figure 6
# load the following objects (see XX for more information on how these were generated)
seurat <- readRDS("Hybrid_seurat_3.Rds")

### Figure 6G
# UMAP plot
DimPlot(seurat)

### Figure 6H
# DotPlot of cluster-specific and osteoblast or adipocyte representative genes
DotPlot(seurat, features=c("MDM2", "SLC7A11", "ZMAT3", "PDGFRA", "PDGFRB", "FGF10", "ADIPOQ", "PLIN1", "HSD11B1", "POSTN", "RUNX2", "COL1A1", "Ob_scaled", "Ad_scaled", "LMCD1", "ALPL", "PPARG"))

### Figure 6I
# ViolinPlot of scaled values for osteoblast genes across cell clusters
VlnPlot(seurat, "Ob_scaled")

# ViolinPlot of scaled values for adipocyte genes across cell clusters
VlnPlot(seurat, "Ad_scaled")

### Figure 6J
# Ob scaled vs. Ad scaled for each cluster together with Doublet score
h <- "0" #insert cluster of interest, 0, 1, 2 or 3
x <- subset(subP, cells=rownames(subP@meta.data[subP@meta.data$seurat_clusters %in% h,]))
df <- cbind(as.data.frame(x$DoubletScore), as.data.frame(x$Ad_scaled), as.data.frame(x$Ob_scaled))
colnames(df) <- c("DoubletScore", "Ad_scaled", "Ob_scaled")

#Plot figure 6J
ggplot() + 
  geom_point(data=df, aes(x=Ad_scaled, y=Ob_scaled, colour=DoubletScore)) + 
  xlim(-60,480) + ylim(-70,150) + 
  scale_colour_gradient(low="lightgrey", high="blue", na.value="red") 
