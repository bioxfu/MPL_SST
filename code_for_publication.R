suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

## Preprocess
seurat_obj <- Load10X_Spatial('/path/of/spaceranger/output') 
seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)

cutoff <- 1.5
x <- seurat_obj[['SCT']]@data@x
seurat_obj[['SCT']]@data@x[x < cutoff] <- 0
seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = FALSE)
seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = 1:30, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, verbose = FALSE, resolution = 0.6)
n <- length(unique(seurat_obj$seurat_clusters))


## Figure 1_I
DimPlot(seurat_obj, reduction = "tsne", label = TRUE, label.color = 'white', 
        label.box = TRUE, repel = TRUE, cols = DiscretePalette(n)) + 
  scale_fill_manual(values = DiscretePalette(n))
SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3, pt.size.factor = 1.3, repel = TRUE) +
  scale_fill_manual(values = DiscretePalette(n))


## Figure 1_J 
markers_c13 <- FindMarkers(seurat_obj, ident.1 = '13', test.use = 'MAST', 
                           features = VariableFeatures(seurat_obj), only.pos = TRUE, 
                           min.diff.pct = 0.2)
markers_c13_top <- markers_c13 %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(!grepl('transgene-', gene)) %>% 
  filter(!grepl('Hb[ab]-', gene)) %>% 
  filter(avg_logFC > 1, p_val_adj < 0.05) %>% 
  top_n(n = 5, wt = avg_logFC)
seurat_obj$target_cluster <- as.character(seurat_obj$seurat_clusters)
seurat_obj$target_cluster[seurat_obj$seurat_clusters==13] <- '13'
seurat_obj$target_cluster[seurat_obj$seurat_clusters!=13] <- 'Other'

Idents(seurat_obj) <- seurat_obj$target_cluster
seurat_obj_mean <- AverageExpression(seurat_obj, assays = 'SCT', slot = 'data', return.seurat = TRUE)[['SCT']]
seurat_obj_mean <- as.data.frame(seurat_obj_mean@data)

plot(seurat_obj_mean$Other, seurat_obj_mean$`13`, cex=1.5, col = 'gray', pch = 19, ylim=c(0,7), 
     xlim = c(0,7), ylab='Average expression in #13 cluster', 
     xlab = 'Average expression in other clusters')
points(seurat_obj_mean[markers_c13_top$gene,], cex=1.5,pch = 19, col = 'red')
text(seurat_obj_mean[markers_c13_top$gene,], markers_c13_top$gene, pos = 4, col = 'red')


## Figure 1_K
SpatialFeaturePlot(seurat_obj, features = c("Sst", "Pth2"), alpha = c(0.1, 1), ncol = 2)


## Figure 1_L
p1 <- FeaturePlot(seurat_obj, features = c("Sst"), reduction = 'tsne') + NoLegend()
p2 <- FeaturePlot(seurat_obj, features = c("Pth2"), reduction = 'tsne')
p1 + p2 + plot_layout(ncol = 2)


## Figure S2_H
VlnPlot(seurat_obj, features = c("Sst", "Pth2", "Slc17a6", 'Slc32a1', 'Gad1', 'Gad2'), 
        pt.size = 0, cols = DiscretePalette(n), ncol = 3, same.y.lims = TRUE)

