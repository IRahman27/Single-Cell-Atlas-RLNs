# Are there batch effects?

library(Seurat)
library(SeuratObject)
library(harmony)

# UMAP by computed clusters before integration
p1 <- DimPlot(RLN, reduction = "umap", group.by = 'seurat_clusters') 

# UMAP by Origin
p2 <- DimPlot(RLN, reduction = "umap", group.by="orig.ident") 

# UMAP by sequencing type (3' or 5')
# Create a vector of sequencing methods based on orig.ident
sequencing_methods <- ifelse(RLN$orig.ident == "Kim", "3'",
                             ifelse(RLN$orig.ident == "Okosun", "5'",
                                    ifelse(RLN$orig.ident == "Roider", "3'", NA)))

# Add the sequencing_method column to metadata
RLNy@meta.data$sequencing_method <- sequencing_methods

p3 <- DimPlot(RLN, reduction = "umap", group.by ="sequencing_method") 

# Look at UMAPs together
p1 + p2 + p3

# UMAP by computed clusters before integration
p4 <- DimPlot(RLN, reduction = "umap", group.by = 'seurat_clusters') 

# UMAP by Origin
p5 <- DimPlot(RLN, reduction = "umap", group.by="orig.ident") 

# UMAP by sequencing type (3' or 5')
# Create a vector of sequencing methods based on orig.ident
sequencing_methods <- ifelse(RLN$orig.ident == "Kim", "3'",
                             ifelse(RLN$orig.ident == "Okosun", "5'",
                                    ifelse(RLN$orig.ident == "Roider", "3'", NA)))

# Add the sequencing_method column to metadata
RLN@meta.data$sequencing_method <- sequencing_methods

p6 <- DimPlot(RLN, reduction = "umap", group.by ="seurat_clusters") 

# Look at UMAPs together
p4 + p5 + p6 + p1 + p2 + p3

grid.arrange(p4, p5, p6, p1, p2, p3, ncol = 3)

# Integration

Integrating the data using Harmony

```{r}

RLN <- RLN %>%
  RunHarmony(group.by.vars = c('orig.ident', 'sequencing_method'), theta = c(2, 6), plot_convergence = FALSE)



RLN.harmony.embed <- Embeddings(RLN, "harmony")
RLN.harmony.embed[1:10,1:10]


# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
RLN <- RLN %>%
  RunUMAP(reduction = 'harmony', dims = 1:10) %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>%
  FindClusters(resolution = 0.5)

# Visualise

p7 <- DimPlot(RLN, reduction = 'umap', label = TRUE)
p8 <- DimPlot(RLN, reduction = 'umap', group.by = 'orig.ident')
p9 <- DimPlot(RLN, reduction = 'umap', group.by = 'sequencing_method')

# Before integration vs after integration using Harmony 
p7 + p8 + p9

RLN <- RLN.Integrated.Harmony