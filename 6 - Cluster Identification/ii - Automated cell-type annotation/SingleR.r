# SingleR

library(SingleR)
library(SingleCellExperiment)
library(celldex)

# Load the relevant reference dataset
ref <- celldex::DatabaseImmuneCellExpressionData()

# Count data from RLN object
RLN.counts <- GetAssayData(RLN, slot = 'counts')

# Run Single R
SR.results <- SingleR(test = RLN.counts, 
                             ref = ref,
                             labels = ref$label.main)


# Add the SingleR labels to RLN metadata
RLN$singleR.labels <- SR.results$labels[match(rownames(RLN@meta.data), rownames(SR.results))]
s1 <- DimPlot(RLN, reduction = 'umap', group.by = 'singleR.labels', label = TRUE)

# Run with another relevant reference data set for comparison
ref <- celldex::HumanPrimaryCellAtlasData()

# Run SingleR with the new reference dataset
SR.results <- SingleR(test = RLN.counts, 
                             ref = ref,
                             labels = ref$label.main)


# Add the SingleR labels for this reference
RLN$singleR.labels <- SR.results$labels[match(rownames(RLN@meta.data), rownames(SR.results))]

s2 <- DimPlot(RLN, reduction = 'umap', group.by = 'singleR.labels', label = TRUE)

# Visualise both 
s2 + s1