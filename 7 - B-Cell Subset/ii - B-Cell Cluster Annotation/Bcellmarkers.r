# Finding markers for the B-Cell subset:

library(Seurat)
library(SeuratObject)

BCells.markers <- FindAllMarkers(BCells)

BCells.markers <- BCells.markers %>%
                      filter(p_val_adj < 0.05, 
                            avg_log2FC > 0.3)

BCells.markers.list <- BCells.markers %>%
                    group_by(cluster) %>%
                    slice_max(n = 200, order_by = avg_log2FC)



write.xlsx(BCells.markers.list, "BCells_top_markers.xlsx")

BCells.markers <- FindAllMarkers(BCells)

BCells.markers <- BCells.markers %>%
                      filter(p_val_adj < 0.05, 
                            avg_log2FC > 0.3)

BCells.markers.list2 <- BCells.markers %>%
                    group_by(cluster) %>%
                    slice_max(n = 5, order_by = avg_log2FC)


write.xlsx(BCells.markers.list, "BCells_top_markers.xlsx")

# Feature plots to identify the specific cell populations in RLN_BCells -  these are marker genes from published litertature

# Cycling B-Cells
FeaturePlot(RLN_BCells2, features = c("HMGB2",	"TUBA1B",	"MKI67", "UBEC2C",	"AURKB", "PCNA", 	"MKI67",	"CDK1",	"CDC20"))

# LZ, DZ
FeaturePlot(RLN_BCells2, features = c("CD83",	"BCL2A1", "CD40", "CXCR4",	"AICDA",	"FOXO1"))

# Naive (2), activated (5)
FeaturePlot(RLN_BCells2, features = c("CD72", "IGHM", "CD69",	"CD83",	"JUN",	"DUSP2",	"EGR1"))

# CD20+ B Cells
FeaturePlot(RLN_BCells2, features = c("MS4A1"))

# Plasma Cells
FeaturePlot(RLN_BCells2, features = c("SDC1", "CD27", "JCHAIN", "IGHG1", "IGHA1", "IGHM", "IGKC", "MZB1"))

# Memory B Cells
FeaturePlot(RLN_BCells2, features = c("MS4A1", "BANK1", "CD27", "CD74", "B2M", "IGHG1", "GPR183",	"CD44",	"KLF2",	"TNFRSF13B"))

#Naive
FeaturePlot(RLN_BCells2, features = c("TCL1A", "IGHD", "IGHM", "IGKC", "CD74", "TXNIP", "BANK1",	"SELL","FCER2", "FCMR", "CD72",	"IGHM"))		

# Class switch recombination genes
FeaturePlot(RLN_BCells2, features = c("APEX1", 	"NME2", "NCL",	"DDX21"))

# Marginal Zone B Cels
FeaturePlot(RLN_BCells2, features = c("CD1C", "CD21", "CD45RB"))

#  FCRL3
FeaturePlot(BCells, features = c("CD19",	"IGHD",	"CD38"))

# Germinal 
FeaturePlot(BCells, features = c("POU2AF1", "CD40", "SUGCT", "CD38", "CD27",	"CD19", "MS4A1", "MME", "CD10"))

#Transitional 
FeaturePlot(BCells, features = c("CD24", "MYO1C", "MS4A1"))

# DZ
FeaturePlot(BCells, features = c("CXCR4", "BCL6", "AICDA", "MYC", "PAX5", "FOX01"))

#LZ
FeaturePlot(BCells, features = c("CD83", "CD86", "CD40", "CXCR4", "AID"))


# Identifying the B-Cell groups 

BCells_final <- BCells
BCells_final.markers <- FindAllMarkers(BCells_final)
BCells_final.markers <- BCells_final.markers %>%
                      filter(p_val_adj < 0.05, 
                            avg_log2FC > 0.3)

# Filtering markers based on p-value and avg log2FC
BCells_final.markers <- BCells_final.markers %>%
  filter(p_val_adj < 0.05, 
         avg_log2FC > 0.3) %>%
  filter(!str_detect(gene, "^MT"))

# Return the rop 15 markers for each cluster
BCells_final.markers.list <- BCells_final.markers %>%
                    group_by(cluster) %>%
                    slice_max(n = 15, order_by = avg_log2FC)

# Visualise this in a heatmap
h1 <- DoHeatmap(BCells_final, features = BCells_final.markers.list$gene) 
```
# Rename the clusters
new.cluster.ids <- c("Naive", "MBC", "MBC", "Activated", "Naive", "MBC", "MBC", "MBC", "MBC", "Plasma Cells", "PreGC", "GC", "GC", "Plasma Cells", "PreGC", "GC")
BCells_final <- BCells
names(new.cluster.ids) <- levels(BCells_final)
BCells_final <- RenameIdents(BCells_final, new.cluster.ids)

