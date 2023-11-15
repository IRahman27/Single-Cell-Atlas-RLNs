

# DOUBLET PIPELINE

# First I am going to filter out all the low quality cells and pre-process each dataset individually. Then I will be able run DoubletFinder and remove the doublets. 
# Once this is done, I will re-run the Seurat pipeline with the doublets removed.

# DoubletFinder requires processing each sample individually through to UMAP:
library(Seurat)
library(Seurat_Object)
library(doubletFinder_v3)

# Kim01
#QC + Filtering
Kim01[["percent.mt"]] <- PercentageFeatureSet(Kim01, pattern = "^MT-")
VlnPlot(Kim01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Kim01 <- subset(Kim01, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
#Pre-processing
# normalise and run dimensionality reduction on control dataset
Kim01 <- SCTransform(Kim01, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

# Kim02
# QC + Filtering
Kim02[["percent.mt"]] <- PercentageFeatureSet(Kim02, pattern = "^MT-")
VlnPlot(Kim02, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Kim02 <- subset(Kim02, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Kim02 <- SCTransform(Kim02, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

# Kim03
# QC + Filtering
Kim03[["percent.mt"]] <- PercentageFeatureSet(Kim03, pattern = "^MT-")
VlnPlot(Kim03, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Kim03 <- subset(Kim03, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Kim03 <- SCTransform(Kim03, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

# Kim04
# QC + Filtering
Kim04[["percent.mt"]] <- PercentageFeatureSet(Kim04, pattern = "^MT-")
VlnPlot(Kim04, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Kim04 <- subset(Kim04, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Kim04 <- SCTransform(Kim04, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

# Kim05
# QC + Filtering
Kim05[["percent.mt"]] <- PercentageFeatureSet(Kim05, pattern = "^MT-")
VlnPlot(Kim05, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Kim05 <- subset(Kim05, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Kim05 <- SCTransform(Kim05, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

# Kim06
# QC + Filtering
Kim06[["percent.mt"]] <- PercentageFeatureSet(Kim06, pattern = "^MT-")
VlnPlot(Kim06, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Kim06 <- subset(Kim06, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Kim06 <- SCTransform(Kim06, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

#Kim07
# QC + Filtering
Kim07[["percent.mt"]] <- PercentageFeatureSet(Kim07, pattern = "^MT-")
VlnPlot(Kim07, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Kim07 <- subset(Kim07, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Kim07 <- SCTransform(Kim07, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

# Kim08
# QC + Filtering
Kim08[["percent.mt"]] <- PercentageFeatureSet(Kim08, pattern = "^MT-")
VlnPlot(Kim08, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Kim08 <- subset(Kim08, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Kim08 <- SCTransform(Kim08, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

#Kim09
# QC + Filtering
Kim09[["percent.mt"]] <- PercentageFeatureSet(Kim09, pattern = "^MT-")
VlnPlot(Kim09, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Kim09 <- subset(Kim09, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Kim09 <- SCTransform(Kim09, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

#Kim10
# QC + Filtering
Kim10[["percent.mt"]] <- PercentageFeatureSet(Kim10, pattern = "^MT-")
VlnPlot(Kim10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Kim10 <- subset(Kim10, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Kim10 <- SCTransform(Kim10, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

#Okosun01
# QC + Filtering
Okosun01[["percent.mt"]] <- PercentageFeatureSet(Okosun01, pattern = "^MT-")
VlnPlot(Okosun01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Okosun01 <- subset(Okosun01, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Okosun01 <- SCTransform(Okosun01, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

#Okosun02
# QC + Filtering
Okosun02[["percent.mt"]] <- PercentageFeatureSet(Okosun02, pattern = "^MT-")
VlnPlot(Okosun02, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Okosun02 <- subset(Okosun02, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Okosun02 <- SCTransform(Okosun02, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

#Okosun03
# QC + Filtering
Okosun03[["percent.mt"]] <- PercentageFeatureSet(Okosun03, pattern = "^MT-")
VlnPlot(Okosun03, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Okosun03 <- subset(Okosun03, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Okosun03 <- SCTransform(Okosun01, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

#Roider01
# QC + Filtering
Roider01[["percent.mt"]] <- PercentageFeatureSet(Roider01, pattern = "^MT-")
VlnPlot(Roider01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Roider01 <- subset(Roider01, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Roider01 <- SCTransform(Roider01, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

#Roider02
# QC + Filtering
Roider02[["percent.mt"]] <- PercentageFeatureSet(Roider02, pattern = "^MT-")
VlnPlot(Roider02, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Roider02 <- subset(Roider02, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Roider02 <- SCTransform(Roider02, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

#Roider03
# QC + Filtering
Roider03[["percent.mt"]] <- PercentageFeatureSet(Roider03, pattern = "^MT-")
VlnPlot(Roider03, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Roider03 <- subset(Roider03, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# Pre-processing
# normalise and run dimensionality reduction on control dataset
Roider03 <- SCTransform(Roider03, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

# Now that each sample has been processed through to UMAP, DoubletFinder can be run on the samples.


# Number of cells in Kim01
length(row.names(Kim01@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 2.4%

# pK Identification (no ground-truth)
sweep.res.list_Kim01 <- paramSweep_v3(Kim01, PCs = 1:20, sct = TRUE)
sweep.stats_Kim01 <- summarizeSweep(sweep.res.list_Kim01, GT = FALSE)
bcmvn_Kim01 <- find.pK(sweep.stats_Kim01)

pK <- bcmvn_Kim01 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  summarize(pK)
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Kim01@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.024*nrow(Kim01@meta.data))  # Assuming 2.4% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Kim01 <- doubletFinder_v3(Kim01, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)


# Number of cells in Kim02
length(row.names(Kim02@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 2.4%

# pK Identification (no ground-truth)
sweep.res.list_Kim02 <- paramSweep_v3(Kim02, PCs = 1:20, sct = TRUE)
sweep.stats_Kim02 <- summarizeSweep(sweep.res.list_Kim02, GT = FALSE)
bcmvn_Kim02 <- find.pK(sweep.stats_Kim02)

pK <- bcmvn_Kim02 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  summarize(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Kim02@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.024*nrow(Kim02@meta.data))  # Assuming 2.4% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Kim02 <- doubletFinder_v3(Kim02, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)




# Number of cells in Kim03
length(row.names(Kim03@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 2.4%

# pK Identification (no ground-truth)
sweep.res.list_Kim03 <- paramSweep_v3(Kim03, PCs = 1:20, sct = TRUE)
sweep.stats_Kim03 <- summarizeSweep(sweep.res.list_Kim03, GT = FALSE)
bcmvn_Kim03 <- find.pK(sweep.stats_Kim03)


pK <- bcmvn_Kim03 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Kim03@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.024*nrow(Kim03@meta.data))  # Assuming 2.4% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# Run DoubletFinder 
Kim03 <- doubletFinder_v3(Kim03, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)



# Number of cells in Kim04
length(row.names(Kim04@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 2.4%

# pK Identification (no ground-truth)
sweep.res.list_Kim04 <- paramSweep_v3(Kim04, PCs = 1:20, sct = TRUE)
sweep.stats_Kim04 <- summarizeSweep(sweep.res.list_Kim04, GT = FALSE)
bcmvn_Kim04 <- find.pK(sweep.stats_Kim04)

pK <- bcmvn_Kim04 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Kim04@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.024*nrow(Kim04@meta.data))  # Assuming 2.4% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Kim04 <- doubletFinder_v3(Kim04, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)


# Number of cells in Kim05
length(row.names(Kim05@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 2.4%

# pK Identification (no ground-truth)
sweep.res.list_Kim05 <- paramSweep_v3(Kim05, PCs = 1:20, sct = TRUE)
sweep.stats_Kim05 <- summarizeSweep(sweep.res.list_Kim05, GT = FALSE)
bcmvn_Kim05 <- find.pK(sweep.stats_Kim05)

pK <- bcmvn_Kim05 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

# Homotypic Doublet Proportion Estimate
annotations <- Kim05@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.024*nrow(Kim05@meta.data))  # Assuming 2.4% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Kim05 <- doubletFinder_v3(Kim05, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)


# Number of cells in Kim06
length(row.names(Kim06@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 2.4%

# pK Identification (no ground-truth)
sweep.res.list_Kim06 <- paramSweep_v3(Kim06, PCs = 1:20, sct = TRUE)
sweep.stats_Kim06 <- summarizeSweep(sweep.res.list_Kim06, GT = FALSE)
bcmvn_Kim06 <- find.pK(sweep.stats_Kim06)

pK <- bcmvn_Kim06 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Kim06@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.024*nrow(Kim06@meta.data))  # Assuming 2.4% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Kim06 <- doubletFinder_v3(Kim06, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)


# Number of cells in Kim07
length(row.names(Kim07@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 4.8%

# pK Identification (no ground-truth)
sweep.res.list_Kim07 <- paramSweep_v3(Kim07, PCs = 1:20, sct = TRUE)
sweep.stats_Kim07 <- summarizeSweep(sweep.res.list_Kim07, GT = FALSE)
bcmvn_Kim07 <- find.pK(sweep.stats_Kim07)

pK <- bcmvn_Kim07 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Kim07@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.048*nrow(Kim07@meta.data))  # Assuming 4.8% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Kim07 <- doubletFinder_v3(Kim07, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)



# Number of cells in Kim08
length(row.names(Kim08@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 3.2%

# pK Identification (no ground-truth)
sweep.res.list_Kim08 <- paramSweep_v3(Kim08, PCs = 1:20, sct = TRUE)
sweep.stats_Kim08 <- summarizeSweep(sweep.res.list_Kim08, GT = FALSE)
bcmvn_Kim08 <- find.pK(sweep.stats_Kim08)

pK <- bcmvn_Kim08 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Kim08@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.032*nrow(Kim08@meta.data))  # Assuming 3.2% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Kim08 <- doubletFinder_v3(Kim08, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)



# Number of cells in Kim09
length(row.names(Kim09@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 3.2%

# pK Identification (no ground-truth)
sweep.res.list_Kim09 <- paramSweep_v3(Kim09, PCs = 1:20, sct = TRUE)
sweep.stats_Kim09 <- summarizeSweep(sweep.res.list_Kim09, GT = FALSE)
bcmvn_Kim09 <- find.pK(sweep.stats_Kim09)


pK <- bcmvn_Kim09 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Kim09@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.032*nrow(Kim09@meta.data))  # Assuming 3.2% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Kim09 <- doubletFinder_v3(Kim09, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)


# Number of cells in Kim10
length(row.names(Kim10@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 4%

# pK Identification (no ground-truth)
sweep.res.list_Kim10 <- paramSweep_v3(Kim10, PCs = 1:20, sct = TRUE)
sweep.stats_Kim10 <- summarizeSweep(sweep.res.list_Kim10, GT = FALSE)
bcmvn_Kim10 <- find.pK(sweep.stats_Kim10)

pK <- bcmvn_Kim10 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  summarize(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Kim10@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.04*nrow(Kim10@meta.data))  # Assuming 4% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Kim10 <- doubletFinder_v3(Kim10, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)

# Number of cells in Okosun01
length(row.names(Okosun01@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 7.2%

# pK Identification (no ground-truth)
sweep.res.list_Okosun01 <- paramSweep_v3(Okosun01, PCs = 1:20, sct = TRUE)
sweep.stats_Okosun01 <- summarizeSweep(sweep.res.list_Okosun01, GT = FALSE)
bcmvn_Okosun01 <- find.pK(sweep.stats_Okosun01)

pK <- bcmvn_Okosun01 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  summarize(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Okosun01@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.072*nrow(Okosun01@meta.data))  # Assuming 7.2% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Okosun01 <- doubletFinder_v3(Okosun01, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)



# Number of cells in Okosun02
length(row.names(Okosun02@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 5.6%

# pK Identification (no ground-truth)
sweep.res.list_Okosun02 <- paramSweep_v3(Okosun02, PCs = 1:20, sct = TRUE)
sweep.stats_Okosun02 <- summarizeSweep(sweep.res.list_Okosun02, GT = FALSE)
bcmvn_Okosun02 <- find.pK(sweep.stats_Okosun02)

pK <- bcmvn_Okosun02 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  summarize(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Okosun02@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.056*nrow(Okosun02@meta.data))  # Assuming 5.6% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Okosun02 <- doubletFinder_v3(Okosun02, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)




# Number of cells in Okosun03
length(row.names(Okosun03@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 4.8%

# pK Identification (no ground-truth)
sweep.res.list_Okosun03 <- paramSweep_v3(Okosun03, PCs = 1:20, sct = TRUE)
sweep.stats_Okosun03 <- summarizeSweep(sweep.res.list_Okosun03, GT = FALSE)
bcmvn_Okosun03 <- find.pK(sweep.stats_Okosun03)

pK <- bcmvn_Okosun03 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  summarize(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Okosun03@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.048*nrow(Okosun03@meta.data))  # Assuming 4.8% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Okosun03 <- doubletFinder_v3(Okosun03, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)


# Number of cells in Roider01
length(row.names(Roider01@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 2.4%

# pK Identification (no ground-truth)
sweep.res.list_Roider01 <- paramSweep_v3(Roider01, PCs = 1:20, sct = TRUE)
sweep.stats_Roider01 <- summarizeSweep(sweep.res.list_Roider01, GT = FALSE)
bcmvn_Roider01 <- find.pK(sweep.stats_Roider01)

pK <- bcmvn_Roider01 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Roider01@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.024*nrow(Roider01@meta.data))  # Assuming 2.4% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Roider01 <- doubletFinder_v3(Roider01, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)


# Number of cells in Roider02
length(row.names(Roider02@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 1.6%

# pK Identification (no ground-truth)
sweep.res.list_Roider02 <- paramSweep_v3(Roider02, PCs = 1:20, sct = TRUE)
sweep.stats_Roider02 <- summarizeSweep(sweep.res.list_Roider02, GT = FALSE)
bcmvn_Roider02 <- find.pK(sweep.stats_Roider02)

pK <- bcmvn_Roider02 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  summarize(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Roider02@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.016*nrow(Roider02@meta.data))  # Assuming 1.6% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Roider02 <- doubletFinder_v3(Roider02, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)


# Number of cells in Roider03
length(row.names(Roider03@meta.data))

# Expected multiplet rate from 10X Genomics website
# ~ 2.4%

# pK Identification (no ground-truth)
sweep.res.list_Roider03 <- paramSweep_v3(Roider03, PCs = 1:20, sct = TRUE)
sweep.stats_Roider03 <- summarizeSweep(sweep.res.list_Roider03, GT = FALSE)
bcmvn_Roider03 <- find.pK(sweep.stats_Roider03)


pK <- bcmvn_Roider03 %>% # select the pK that corresponds to max bcmvn to optimise doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


# Homotypic Doublet Proportion Estimate
annotations <- Roider03@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.024*nrow(Roider03@meta.data))  # Assuming 2.4% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Roider03 <- doubletFinder_v3(Roider03, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = TRUE)


# Renaming the Doublet Columns
colnames(Kim01@meta.data)[colnames(Kim01@meta.data) == "old_column_name"] <- "df_predictions"
colnames(Kim02@meta.data)[colnames(Kim02@meta.data) == "old_column_name"] <- "df_predictions"

