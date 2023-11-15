# Merge and pre-process as a singular object

library(Seurat)
library(SeuratObject)

# Merge into 1 object
RLN <- merge(Kim01, y = c(Kim02, Kim03, Kim04, Kim05, Kim06, Kim07, Kim08, Kim09, Kim10, Okosun01, Okosun02, Okosun03, Roider01, Roider02, Roider03),
             add.cell.ids = c("K1", "K2", "K3", "K4", "K5", "K6", "K7", "K8", "K9", "K10", "O1", "O2", "O3", "R1", "R2", "R3"), project = 'RLNs')

# How many cells from each study after removing doublets?
table(RLN$orig.ident)


# Pre-Processing the merged data
RLN[["percent.mt"]] <- PercentageFeatureSet(RLN, pattern = "^MT-")
VlnPlot(RLN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RLN <- SCTransform(RLN, vst.flavor = "v2", verbose = FALSE) 
RLN <- RunPCA(RLN, dims = 1:20, verbose = FALSE) 
RLN_heatmap <- DimHeatmap(RLN, dims = 1:20,balanced = TRUE)
ElbowPlot(RLN) 

RLN <- FindNeighbors(RLN, reduction = "pca", dims = 1:15, verbose = FALSE) 
RLN <- FindClusters(RLN, resolution = 0.3, verbose = FALSE) 
RLN <- RunUMAP(RLN, reduction = "pca", dims = 1:15, verbose = FALSE) 

# DoubletFinder had different column names for predictions for each object


RLN@meta.data$df_predictions <- c (Kim01@meta.data$DF.classifications_0.25_0.05_76, Kim02@meta.data$DF.classifications_0.25_0.06_66, Kim03@meta.data$DF.classifications_0.25_0.06_64, Kim04@meta.data$DF.classifications_0.25_0.19_62, Kim05@meta.data$DF.classifications_0.25_0.06_67, Kim06@meta.data$DF.classifications_0.25_0.12_59, Kim07@meta.data$DF.classifications_0.25_0.06_242, Kim08@meta.data$DF.classifications_0.25_0.04_124, Kim09@meta.data$DF.classifications_0.25_0.04_107, Kim10@meta.data$DF.classifications_0.25_0.23_158, Okosun01@meta.data$DF.classifications_0.25_0.18_605, Okosun02@meta.data$DF.classifications_0.25_0.11_364, Okosun03@meta.data$DF.classifications_0.25_0.18_403, Roider01@meta.data$DF.classifications_0.25_0.2_62, Roider02@meta.data$DF.classifications_0.25_0.05_28, Roider03@meta.data$DF.classifications_0.25_0.23_60)


# Visualise the doublets

DimPlot(RLN, reduction = 'umap', group.by = 'df_predictions')
DimPlot(RLN, reduction = 'umap', group.by = 'scdoublet_predictions')

# Identifying the doublets

# Access the metadata columns
doublet_finder_col <- "df_predictions"
scrublet_col <- "scdoublet_predictions"

# Get the values from metadata
doublet_finder_values <- RLN@meta.data[[doublet_finder_col]]
scrublet_values <- RLN@meta.data[[scrublet_col]]

# Create a variable for cells identified as doublets by both methods
both_doublets <- doublet_finder_values == "Doublet" & scrublet_values == "Doublet"

# Create a variable for doublets only in DF.classifications_0.25_0.07_71
df_doublets <- doublet_finder_values == "Doublet" & scrublet_values == "Singlet"

# Create a variable for doublets only in scdoublet_predictions
scrublet_doublets <- doublet_finder_values == "Singlet" & scrublet_values == "Doublet"

# Calculate the counts
both_doublets_count <- sum(both_doublets)
df_doublets_count <- sum(df_doublets)
scrublet_doublets_count <- sum(scrublet_doublets)

# Print the results
cat("Cells identified as doublets by both methods:", both_doublets_count, "\n")
cat("Cells identified as doublets only by DoubletFinder:", df_doublets_count, "\n")
cat("Cells identified as doublets only by Scrublet:", scrublet_doublets_count, "\n")

RLN_both <- subset(RLN, subset = df_predictions == "Doublet" & scdoublet_predictions == "Doublet")
RLN_df <- subset(RLN, subset = df_predictions == "Doublet")
RLN_sc <- subset(RLN, scdoublet_predictions == "Doublet")

# Number of cells identified by both tools as doublets
length(rownames(RLN_both@meta.data)) # 736

# Number of cells identified by DF as doublets
length(rownames(RLN_df@meta.data)) # 2547

# Number of cells identified by Scrublet as doublets
length(rownames(RLN_sc@meta.data)) # 1058

# Access gene expression data and gene names
gene_expression <- LayerData(object = RLN_both, assay = "RNA")
gene_names <- rownames(gene_expression)
```
# Biology of the doublets 

# Count matrix for gene expression values
RLN_cm <- LayerData(RLN)

# Top genes for first cell identified as doublet by scrublet
sort.int(RLN_cm[,'K1_ACTGAACCAATGAAAC_LN_01'], decreasing = TRUE)

# Remove the doublets 

# Re-run the seurat pipeline from scratch

# Creating Seurat objects for Kim data
Kim01 <- CreateSeuratObject(Kim01.data, project = 'Kim')
Kim02 <- CreateSeuratObject(Kim02.data, project = 'Kim')
Kim03 <- CreateSeuratObject(Kim03.data, project = 'Kim')
Kim04 <- CreateSeuratObject(Kim04.data, project = 'Kim')
Kim05 <- CreateSeuratObject(Kim05.data, project = 'Kim')
Kim06 <- CreateSeuratObject(Kim06.data, project = 'Kim')
Kim07 <- CreateSeuratObject(Kim07.data, project = 'Kim')
Kim08 <- CreateSeuratObject(Kim08.data, project = 'Kim')
Kim09 <- CreateSeuratObject(Kim09.data, project = 'Kim')
Kim10 <- CreateSeuratObject(Kim10.data, project = 'Kim')

# Create Seurat objects for Okosun data
Okosun01 <- CreateSeuratObject(Okosun01.data, project = 'Okosun')
Okosun02 <- CreateSeuratObject(Okosun02.data, project = 'Okosun')
Okosun03 <- CreateSeuratObject(Okosun03.data, project = 'Okosun')

# Create Seurat objects for Roider data
Roider01 <- CreateSeuratObject(Roider01.data, project = 'Roider')
Roider02 <- CreateSeuratObject(Roider02.data, project = 'Roider')
Roider03 <- CreateSeuratObject(Roider03.data, project = 'Roider')

# Merge into 1 object
RLN <- merge(Kim01, y = c(Kim02, Kim03, Kim04, Kim05, Kim06, Kim07, Kim08, Kim09, Kim10, Okosun01, Okosun02, Okosun03, Roider01, Roider02, Roider03),
             add.cell.ids = c("K1", "K2", "K3", "K4", "K5", "K6", "K7", "K8", "K9", "K10", "O1", "O2", "O3", "R1", "R2", "R3"), project = 'RLNs')


# Identify the cells to be removed
cells_to_remove <- rownames(RLN_sc@meta.data)

# Subset the merged object to remove the doublets
RLN <- subset(RLN, cells = cells_to_remove, invert = TRUE)

# Pre-processing to normalisation

RLN[["percent.mt"]] <- PercentageFeatureSet(RLN, pattern = "^MT-")
VlnPlot(RLN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RLN <- subset(RLN, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
RLN <- SCTransform(RLN, vst.flavor = "v2", assay = 'RNA', new.assay.name = 'SCT', verbose = FALSE)

# Cell Cycle Scoring

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

s.genes <- s.genes[!s.genes %in% "MLF1IP"]
g2m.genes <- g2m.genes[!g2m.genes %in% "MLF1IP"]

RLN <- CellCycleScoring(RLN, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident = TRUE)

# Run SCTransform and regress out cell cycle and mitochondrial genes

RLN <- SCTransform(RLN, vst.flavor = "v2",  assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), verbose = FALSE) 


# Adjust Variable Features
# Remove unwanted genes 

# - IGHV|IGKV|IGLV … removes B cell receptor (BCR) immunogloublin (Ig) heavy chain (H) and light chain (K/L) variable genes from HVGs

# - TRAV|TRBV|TRDV|TRGV … removes T cell receptor (TCR) variable genes (for alpha/beta/delta/gamma TCRs) from HVGs

# - MT-|MTRN … So cells don’t cluster by mitochondrial content (not usually biologically relevant)


# What are my top 10 highly variable features initially?
RLN@assays$SCT

var.feat <- VariableFeatures(RLN)
new.var.feat <- var.feat[grep("^IGHV|^IGKV|^IGLV|^TRAV|^TRBV|^TRDV|^TRGV|^MT-|^MTRN", var.feat, invert = TRUE)]
VariableFeatures(RLN) <- new.var.feat

# Rest of the pre-processing
RLN <- RunPCA(RLN, dims = 1:20, verbose = FALSE) 
ElbowPlot(RLN) 

RLN <- FindNeighbors(RLN, reduction = "pca", dims = 1:10, verbose = FALSE) 
RLN <- FindClusters(RLN, resolution = 0.5, verbose = FALSE) 
RLN <- RunUMAP(RLN, reduction = "pca", dims = 1:10, verbose = FALSE) 


