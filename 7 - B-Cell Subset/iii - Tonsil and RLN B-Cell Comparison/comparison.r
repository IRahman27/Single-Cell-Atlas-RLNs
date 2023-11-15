library(Seurat)
library(SeuratObject)

# Comparing Tonsil B-Cells and RLN B-Cells
# Merge into 1 object for comparison of cells - do the same identfied cell types cluster together?

# Pre-processing
BCells_merged <- merge(BCells2, y = c(tonsil2),
             add.cell.ids = c("RLN", "Tonsil"), project = 'RLNs')

BCells_merged <- FindVariableFeatures(BCells_merged)
BCells_merged <- ScaleData(BCells_merged)

BCells_merged <- RunPCA(BCells_merged, dims = 1:20, verbose = FALSE) 
ElbowPlot(BCells_merged) 

BCells_merged <- FindNeighbors(BCells_merged, reduction = "pca", dims = 1:10, verbose = FALSE) 
BCells_merged <- FindClusters(BCells_merged, resolution = 0.5, verbose = FALSE) 
BCells_merged <- RunUMAP(BCells_merged, reduction = "pca", dims = 1:10, verbose = FALSE) 

# Visualise merged object
DimPlot(BCells_merged, reduction = "umap", label = TRUE)

#Rename tonsil levels
tonsil2 <- tonsil
new.cluster.ids <- c("T_Naive", "T_Activated", "T_preGC", "T_LZ GC", "T_GC", "T_DZ GC", "T_FCRL2/3high GC", "T_prePB", "T_Plasmablast",   "T_MBC", "T_MBC FCRL4+", "T_Cycling")
names(new.cluster.ids) <- levels(tonsil2)
tonsil2<- RenameIdents(tonsil2, new.cluster.ids)

#Rename RLN levels
BCells2 <- BCells_final
new.cluster.ids <- c("RLN_Naive", "RLN_csMBC", "RLN_Activated", "RLN_ncsMBC", "RLN_Plasma Cells", "RLN_PreGC", "RLN_GC" )
names(new.cluster.ids) <- levels(BCells2)
BCells2<- RenameIdents(BCells2, new.cluster.ids)

# Visualise B-Cell type proportions for RLN and Tonsil

library(ggplot2)
library(dplyr)

# Calculate the frequencies for BCells_final and tonsil
bcells_freq <- table(Idents(BCells_final)) / ncol(BCells_final)*100
tonsil_freq <- table(Idents(tonsil_new)) / ncol(tonsil_new)*100

# Create data frames for frequencies and Seurat object names
bcells_df <- data.frame(Cell_Type = names(bcells_freq), Frequency = bcells_freq, Seurat_Object = "RLN")
tonsil_df <- data.frame(Cell_Type = names(tonsil_freq), Frequency = tonsil_freq, Seurat_Object = "Tonsil")

# Combine the data frames
combined_counts <- rbind(bcells_df, tonsil_df)

# Convert Frequency to a numeric value
combined_counts$Frequency.Freq <- as.numeric(combined_counts$Frequency.Freq)

# Create the frequency plot
frequency_plot <- ggplot(combined_counts, aes(x = (combined_counts$Frequency.Freq)*100, y = reorder(Cell_Type, combined_counts$Frequency.Freq), fill = Cell_Type)) +
  geom_col(show.legend = FALSE) +
  facet_grid(Seurat_Object~., scales = "free_y", space = "free_y") +
  scale_x_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 60), expand = c(0, 0)) +
  labs(title = "Cell Type Proportions Comparison",
       x = "Percentage",
       y = "Cell Type",
       fill = "Cell Type") +
  theme_minimal()

print(frequency_plot)
