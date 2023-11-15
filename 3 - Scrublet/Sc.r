
# SCRUBLET PIPELINE

# Reticulate so I can use Python
library(reticulate)
use_condaenv("r-reticulate", conda = "auto")

# Installing relevant packages
py_run_string("import pip")
py_run_string("pip.main(['install', 'pandas'])")
py_run_string("pip.main(['install', 'scanpy'])")
py_run_string("pip.main(['install', 'scrublet'])")
py_run_string("pip.main(['install', 'numpy'])")

# Importing the necessary packages 
pd <- import('pandas')
sp <- import('scanpy')
sc <- import('scrublet')
np <- import('numpy')

# Count matrix
Kim01.cm <-LayerData(Kim01, assay='RNA')

# Transposing so cells are rows and genes are columnds
Kim01.cm <- t(as.data.frame(Kim01.cm))

# Initialising scrublet object
K1scrub = sc$Scrublet(Kim01.cm, expected_doublet_rate = 0.024)

Kim01.sc <- K1scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Kim02.cm <-LayerData(Kim02, assay='RNA')

# Transposing so cells are rows and genes are columnds
Kim02.cm <- t(as.data.frame(Kim02.cm))

# Initialising scrublet object
K2scrub = sc$Scrublet(Kim02.cm, expected_doublet_rate = 0.024)

Kim02.sc <- K2scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Kim03.cm <-LayerData(Kim03, assay='RNA')

# Transposing so cells are rows and genes are columnds
Kim03.cm <- t(as.data.frame(Kim03.cm))

# Initialising scrublet object
K3scrub = sc$Scrublet(Kim03.cm, expected_doublet_rate = 0.024)

Kim03.sc <- K3scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Kim04.cm <-LayerData(Kim04, assay='RNA')

# Transposing so cells are rows and genes are columnds
Kim04.cm <- t(as.data.frame(Kim04.cm))

# Initialising scrublet object
K4scrub = sc$Scrublet(Kim04.cm, expected_doublet_rate = 0.024)

Kim04.sc <- K4scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Kim05.cm <-LayerData(Kim05, assay='RNA')

# Transposing so cells are rows and genes are columnds
Kim05.cm <- t(as.data.frame(Kim05.cm))

# Initialising scrublet object
K5scrub = sc$Scrublet(Kim05.cm, expected_doublet_rate = 0.024)

Kim05.sc <- K5scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Kim06.cm <-LayerData(Kim06, assay='RNA')

# Transposing so cells are rows and genes are columnds
Kim06.cm <- t(as.data.frame(Kim06.cm))

# Initialising scrublet object
K6scrub = sc$Scrublet(Kim06.cm, expected_doublet_rate = 0.024)

Kim06.sc <- K6scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Kim07.cm <-LayerData(Kim07, assay='RNA')

# Transposing so cells are rows and genes are columnds
Kim07.cm <- t(as.data.frame(Kim07.cm))

# Initialising scrublet object
K7scrub = sc$Scrublet(Kim07.cm, expected_doublet_rate = 0.048)

Kim07.sc <- K7scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Kim08.cm <-LayerData(Kim08, assay='RNA')

# Transposing so cells are rows and genes are columnds
Kim08.cm <- t(as.data.frame(Kim08.cm))

# Initialising scrublet object
K8scrub = sc$Scrublet(Kim08.cm, expected_doublet_rate = 0.032)

Kim08.sc <- K8scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Kim09.cm <-LayerData(Kim09, assay='RNA')

# Transposing so cells are rows and genes are columnds
Kim09.cm <- t(as.data.frame(Kim09.cm))

# Initialising scrublet object
K9scrub = sc$Scrublet(Kim09.cm, expected_doublet_rate = 0.032)

Kim09.sc <- K9scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Kim10.cm <-LayerData(Kim10, assay='RNA')

# Transposing so cells are rows and genes are columnds
Kim10.cm <- t(as.data.frame(Kim10.cm))

# Initialising scrublet object
K10scrub = sc$Scrublet(Kim10.cm, expected_doublet_rate = 0.04)

Kim10.sc <- K10scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Okosun01.cm <-LayerData(Okosun01, assay='RNA')

# Transposing so cells are rows and genes are columnds
Okosun01.cm <- t(as.data.frame(Okosun01.cm))

# Initialising scrublet object
O1scrub = sc$Scrublet(Okosun01.cm, expected_doublet_rate = 0.072)

Okosun01.sc <- O1scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Okosun02.cm <-LayerData(Okosun02, assay='RNA')

# Transposing so cells are rows and genes are columnds
Okosun02.cm <- t(as.data.frame(Okosun02.cm))

# Initialising scrublet object
O2scrub = sc$Scrublet(Okosun02.cm, expected_doublet_rate = 0.056)

Okosun02.sc <- O2scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Okosun03.cm <-LayerData(Okosun03, assay='RNA')

# Transposing so cells are rows and genes are columnds
Okosun03.cm <- t(as.data.frame(Okosun03.cm))

# Initialising scrublet object
O3scrub = sc$Scrublet(Okosun03.cm, expected_doublet_rate = 0.048)

Okosun03.sc <- O3scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Roider01.cm <-LayerData(Roider01, assay='RNA')

# Transposing so cells are rows and genes are columnds
Roider01.cm <- t(as.data.frame(Roider01.cm))

# Initialising scrublet object
R1scrub = sc$Scrublet(Roider01.cm, expected_doublet_rate = 0.024)

Roider01.sc <- R1scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Count matrix
Roider02.cm <-LayerData(Roider02, assay='RNA')

# Transposing so cells are rows and genes are columnds
Roider02.cm <- t(as.data.frame(Roider02.cm))

# Initialising scrublet object
R2scrub = sc$Scrublet(Roider02.cm, expected_doublet_rate = 0.016)

Roider02.sc <- R2scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

#Roider03
# QC + Filtering
Roider03.s <- Roider03
Roider03.s[["percent.mt"]] <- PercentageFeatureSet(Roider03.s, pattern = "^MT-")
Roider03.s <- subset(Roider03.s, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

# Count matrix
Roider03.cm <-LayerData(Roider03.s, assay='RNA')

# Transposing so cells are rows and genes are columnds
Roider03.cm <- t(as.data.frame(Roider03.cm))

# Initialising scrublet object
R3scrub = sc$Scrublet(Roider03.cm, expected_doublet_rate = 0.024)

Roider03.sc <- R3scrub$scrub_doublets(min_counts=2, 
                                  min_cells=3, 
                                  n_prin_comps= as.integer(30))

# Pre-processing
# normalise and run dimensionality reduction on control dataset
Roider03 <- SCTransform(Roider03, vst.flavor = "v2", verbose = FALSE) %>%
RunPCA(npcs = 20, verbose = FALSE) %>%
RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
FindClusters(resolution = 0.7, verbose = FALSE)

# Adding Scrublet Predictions to the metadata 

Kim01@meta.data$scdoublet_scores <- Kim01.sc[[1]]
Kim01@meta.data$scdoublet_predictions <- Kim01.sc[[2]]

# Rename for uniformity
Kim01$scdoublet_predictions <- ifelse(Kim01$scdoublet_predictions, "Doublet", "Singlet")

Kim02@meta.data$scdoublet_scores <- Kim02.sc[[1]]
Kim02@meta.data$scdoublet_predictions <- Kim02.sc[[2]]

# Rename for uniformity
Kim02$scdoublet_predictions <- ifelse(Kim02$scdoublet_predictions, "Doublet", "Singlet")

Kim03@meta.data$scdoublet_scores <- Kim03.sc[[1]]
Kim03@meta.data$scdoublet_predictions <- Kim03.sc[[2]]

# Rename for uniformity
Kim03$scdoublet_predictions <- ifelse(Kim03$scdoublet_predictions, "Doublet", "Singlet")

Kim04@meta.data$scdoublet_scores <- Kim04.sc[[1]]
Kim04@meta.data$scdoublet_predictions <- Kim04.sc[[2]]

# Rename for uniformity
Kim04$scdoublet_predictions <- ifelse(Kim04$scdoublet_predictions, "Doublet", "Singlet")

Kim05@meta.data$scdoublet_scores <- Kim05.sc[[1]]
Kim05@meta.data$scdoublet_predictions <- Kim05.sc[[2]]

# Rename for uniformity
Kim05$scdoublet_predictions <- ifelse(Kim05$scdoublet_predictions, "Doublet", "Singlet")

Kim06@meta.data$scdoublet_scores <- Kim06.sc[[1]]
Kim06@meta.data$scdoublet_predictions <- Kim06.sc[[2]]

# Rename for uniformity
Kim06$scdoublet_predictions <- ifelse(Kim06$scdoublet_predictions, "Doublet", "Singlet")

Kim07@meta.data$scdoublet_scores <- Kim07.sc[[1]]
Kim07@meta.data$scdoublet_predictions <- Kim07.sc[[2]]

# Rename for uniformity
Kim07$scdoublet_predictions <- ifelse(Kim07$scdoublet_predictions, "Doublet", "Singlet")

Kim08@meta.data$scdoublet_scores <- Kim08.sc[[1]]
Kim08@meta.data$scdoublet_predictions <- Kim08.sc[[2]]

# Rename for uniformity
Kim08$scdoublet_predictions <- ifelse(Kim08$scdoublet_predictions, "Doublet", "Singlet")

Kim09@meta.data$scdoublet_scores <- Kim09.sc[[1]]
Kim09@meta.data$scdoublet_predictions <- Kim09.sc[[2]]

# Rename for uniformity
Kim09$scdoublet_predictions <- ifelse(Kim09$scdoublet_predictions, "Doublet", "Singlet")

Kim10@meta.data$scdoublet_scores <- Kim10.sc[[1]]
Kim10@meta.data$scdoublet_predictions <- Kim10.sc[[2]]

# Rename for uniformity
Kim10$scdoublet_predictions <- ifelse(Kim10$scdoublet_predictions, "Doublet", "Singlet")

Okosun01@meta.data$scdoublet_scores <- Okosun01.sc[[1]]
Okosun01@meta.data$scdoublet_predictions <- Okosun01.sc[[2]]

# Rename for uniformity
Okosun01$scdoublet_predictions <- ifelse(Okosun01$scdoublet_predictions, "Doublet", "Singlet")

Okosun02@meta.data$scdoublet_scores <- Okosun02.sc[[1]]
Okosun02@meta.data$scdoublet_predictions <- Okosun02.sc[[2]]

# Rename for uniformity
Okosun02$scdoublet_predictions <- ifelse(Okosun02$scdoublet_predictions, "Doublet", "Singlet")

Okosun03@meta.data$scdoublet_scores <- Okosun03.sc[[1]]
Okosun03@meta.data$scdoublet_predictions <- Okosun03.sc[[2]]

# Rename for uniformity
Okosun03$scdoublet_predictions <- ifelse(Okosun03$scdoublet_predictions, "Doublet", "Singlet")

Roider01@meta.data$scdoublet_scores <- Roider01.sc[[1]]
Roider01@meta.data$scdoublet_predictions <- Roider01.sc[[2]]

# Rename for uniformity
Roider01$scdoublet_predictions <- ifelse(Roider01$scdoublet_predictions, "Doublet", "Singlet")

Roider02@meta.data$scdoublet_scores <- Roider02.sc[[1]]
Roider02@meta.data$scdoublet_predictions <- Roider02.sc[[2]]

# Rename for uniformity
Roider02$scdoublet_predictions <- ifelse(Roider02$scdoublet_predictions, "Doublet", "Singlet")

Roider03@meta.data$scdoublet_scores <- Roider03.sc[[1]]
Roider03@meta.data$scdoublet_predictions <- Roider03.sc[[2]]

# Rename for uniformity
Roider03$scdoublet_predictions <- ifelse(Roider03$scdoublet_predictions, "Doublet", "Singlet")


