# Single Cell Atlas of Adult Human Reactive Lymph Nodes
Code repository for my postgraduate research project

## Abstract
This project addresses a significant knowledge gap concerning follicular lymphoma (FL) references by focusing on the construction of a single-cell atlas derived from adult human reactive lymph nodes (RLNs) as a potential reference for FL research. The current reference for FL studies is based on paediatric tonsils; however, this study hypothesised that an atlas derived from adult RLNs would offer a more suitable reference due to the immune cell composition's intrinsic characteristics.

The study successfully developed a comprehensive single-cell atlas of adult RLNs, encompassing more than 65,000 individual cells. This atlas was further refined by isolating B-cells, enabling a direct comparison with the B-cell atlas generated from paediatric tonsils. The results demonstrated compelling alignments between the RLN atlas and FL, both in terms of immune cell proportions and gene expression profiles. Notably, the RLN atlas exhibited higher proportions of memory B cells and distinctive gene expression patterns, which closely resembled the characteristics of FL.

The findings underscore the significance of the RLN atlas as a pertinent 'normal' reference for FL investigations. This alignment provides a solid foundation for advancing research on immune cell diversity, interactions, and functional roles within the context of both normal immune responses and disease pathogenesis. However, the study also highlights the challenges associated with constructing a robust single-cell atlas, primarily stemming from the absence of standardized analytical approaches.

The atlas offers opportunities to refine B-cell subtypes, delve into T-cell subsets, and expand its scope by incorporating additional RLN samples. These endeavours would enhance the atlas's depth and applicability. In conclusion, this study has the potential to contribute to our understanding of immune cell dynamics and their roles in health and disease, offering a promising avenue for future research in the realm of immune-related disorders, such as follicular lymphoma.

## Code Implementation
Folders and Content: 

`1 - Data Loading`: Loading data from each study and creating Seurat Objects for each sample

`2 - DoubletFinder`: Running the DoubletFinder pipeline on each sample

`3 - Scrublet`: Running the Scrublet pipeline on each sample

`4 - QC and Processing`: Merging samples and pre-processing as a single object using the Seurat pipeline

`5 - Integration`: Integrating data from the different samples using Harmony

`6- Cluster Identification`: Manual cell annotation and automated cell annotation using SingleR

`7 - B-Cell Subset`: Subsetting RLN object to identify specific B-Cell populations and compare to existing Tonsil atlas


