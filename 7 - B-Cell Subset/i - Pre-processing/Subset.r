library(Seurat)
library(SeuratObject)

# Subset for B-Cells
BCells <- subset(RLN, idents = c('B-Cells', 'Memory B-Cells', 'Plasma Cells'))

# Re-run pre-processing from normalisation for the B-Cell subset
BCells <- SCTransform(BCells, vst.flavor = "v2",  assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), verbose = FALSE)
BCells <- RunPCA(BCells, dims = 1:20, verbose = FALSE)
ElbowPlot(BCells)
BCells <- FindNeighbors(BCells, reduction = "pca", dims = 1:7, verbose = FALSE)
BCells <- FindClusters(BCells, resolution = 0.5, verbose = FALSE)
BCells <- RunUMAP(BCells, reduction = "pca", dims = 1:7, verbose = FALSE)
 
# Visualise B-Cell subset
DimPlot(BCells, reduction = 'umap', label = TRUE)

