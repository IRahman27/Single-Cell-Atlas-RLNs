# The first step is to seperate all the data into samples per Reactive Lymph Node. 
#Load the Kim, Okosun and Roider data individually:

setwd('/data/BCI-OkosunLab/Ilhan/ProjectData')

#Load Kim data
Kim.data <-readRDS("Kim/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")

# The format for the Kim data is one large data frame with all of the samples from that paper. Therefore, the RLN samples need t be extracted using their sample ID's
# Filter rows based on row name suffix

pattern <- ".*LN_01"  
filtered_columns <- grepl(pattern, colnames(Kim.data))  
Kim01.data <- Kim.data[,filtered_columns]
Kim01.data <- Matrix(as.matrix(Kim01.data), sparse = TRUE)

pattern <- ".*LN_02"  
filtered_columns <- grepl(pattern, colnames(Kim.data))  
Kim02.data <- Kim.data[,filtered_columns]
Kim02.data <- Matrix(as.matrix(Kim02.data), sparse = TRUE)

pattern <- ".*LN_03"  
filtered_columns <- grepl(pattern, colnames(Kim.data))  
Kim03.data <- Kim.data[,filtered_columns]
Kim03.data <- Matrix(as.matrix(Kim03.data), sparse = TRUE)

pattern <- ".*LN_04"  
filtered_columns <- grepl(pattern, colnames(Kim.data))  
Kim04.data <- Kim.data[,filtered_columns]
Kim04.data <- Matrix(as.matrix(Kim04.data), sparse = TRUE)


pattern <- ".*LN_05"  
filtered_columns <- grepl(pattern, colnames(Kim.data))  
Kim05.data <- Kim.data[,filtered_columns]
Kim05.data <- Matrix(as.matrix(Kim05.data), sparse = TRUE)

pattern <- ".*LN_06"  
filtered_columns <- grepl(pattern, colnames(Kim.data))  
Kim06.data <- Kim.data[,filtered_columns]
Kim06.data <- Matrix(as.matrix(Kim06.data), sparse = TRUE)


pattern <- ".*LN_07"  
filtered_columns <- grepl(pattern, colnames(Kim.data))  
Kim07.data <- Kim.data[,filtered_columns]
Kim07.data <- Matrix(as.matrix(Kim07.data), sparse = TRUE)


pattern <- ".*LN_08"  
filtered_columns <- grepl(pattern, colnames(Kim.data))  
Kim08.data <- Kim.data[,filtered_columns]
Kim08.data <- Matrix(as.matrix(Kim08.data), sparse = TRUE)

pattern <- ".*LN_11"  
filtered_columns <- grepl(pattern, colnames(Kim.data))  
Kim09.data <- Kim.data[,filtered_columns]
Kim09.data <- Matrix(as.matrix(Kim09.data), sparse = TRUE)

pattern <- ".*LN_12"  
filtered_columns <- grepl(pattern, colnames(Kim.data))  
Kim10.data <- Kim.data[,filtered_columns]
Kim10.data <- Matrix(as.matrix(Kim10.data), sparse = TRUE)

# Load the Okosun data
Okosun01.data <- Read10X(data.dir = "Okosun/R33651")
Okosun02.data <- Read10X(data.dir = "Okosun/R34021")
Okosun03.data <- Read10X(data.dir = "Okosun/R34382")

# Load the Roider data
Roider01.data <- Read10X(data.dir = "Roider/rLN1")
Roider02.data <- Read10X(data.dir = "Roider/rLN2")
Roider03.data <- Read10X(data.dir = "Roider/rLN3")
```


# Create Seurat Objects for the individual datasets so they're ready for QC:

library(Seurat)
library(SeuratObject)

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


