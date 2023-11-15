# Identifying the cell types of each cluster - Feature plots of main immune cells

library(Seurat)
library(SeuratObject)

# B cells, T cells, NK cells
FeaturePlot(RLN.Integrated.Harmony, features = c("MS4A1", "CD19", "CD79A", "CD79B", "CD3", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "KLRC1", "KLRC2" ))

# Dendritic cells
FeaturePlot(RLN.Integrated.Harmony, features = c("CD1C", "FCER1A", "CLEC10A"))

# Macrophages
FeaturePlot(RLN.Integrated.Harmony, features = c("C1QA", "CD68", "TREM2"))

# Monocytes
FeaturePlot(RLN.Integrated.Harmony, features = c("S100A9", "LYZ", "FCN1"))

# Fibroblasts
FeaturePlot(RLN.Integrated.Harmony, features = c("CCDC80", "COL14A1", "FBLN1", "LTBP4", "MFAP4", "OLFML1", "PCOLCE", "SSC5D", "TIMP2"))


# Plasma cells
FeaturePlot(RLN.Integrated.Harmony, features = c("CD38", "CD138", "CD27", "CD45RA"))
FeaturePlot(RLN.Integrated.Harmony, Idents(11), features =c("MS4A1", "CD3D"))

Cluster_11 <- subset(RLN.Integrated.Harmony, idents = '11')

FeaturePlot(Cluster_11, features =c("MS4A1", "CD3D"))

Find markers for each cluster

RLN.markers <- FindAllMarkers(RLN)

RLN.markers.list <- RLN.markers %>%
                    group_by(cluster) %>%
                    slice_max(n = 15, order_by = avg_log2FC)


write.xlsx(RLN.markers.list, "RLN_top_markers.xlsx")
