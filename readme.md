# Single-cell RNA-seq Data Processing Pipeline

**Associated Publication:** Single-Cell Transcriptomic Analysis of Different Liver Fibrosis Models: Elucidating Molecular Distinctions and Commonalities*  
**Data Sources:**  
- [GSE171904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171904) (Liver fibrosis models: oil/ccl4/bdl treatments)  
- [GSE199638](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199638) (Knockout models: CTRL/KO_MMT)  
- [GSE221481](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221481) (Drug interventions: Thioacetamide/Celecoxib)  

**Author:** Guofei Deng 
**Last Updated:** 2025-07-09  

---

## 1. Pipeline Overview
This repository contains the R script for processing and integrating multiple single-cell RNA-seq datasets from GEO. Key steps:

```mermaid
graph TD
    A[Raw GEO Files] --> B[Reorganize to 10X Format]
    B --> C[Create Seurat Objects]
    C --> D[Calculate QC Metrics]
    D --> E[Filter Cells]
    E --> F[Normalize Data]
    F --> G[Run PCA]

## 2. Quick Start
Prerequisites
R (≥4.2)

Conda (recommended)

Installation

bash
conda create -n sc_analysis r-base=4.2 r-seurat=5.0
conda activate sc_analysis
Rscript -e "install.packages(c('tidyverse','clustree','cowplot','data.table'))"

Run Pipeline
1.Place raw data in these directories:
bash
mkdir -p GSE171904_RAW GSE199638_RAW GSE221481_RAW
# Download *.gz files from GEO into corresponding folders

2.Execute the script:
bash
Rscript 1.merge.R

## 3. Output Structure
text
results/
├── GSE171904.rds               # Integrated Seurat object
├── GSE199638.rds
├── GSE221481.rds
└── figures/
    ├── 1.Vlnplot1.pdf          # QC metrics
    ├── 2.Vlnplot2.pdf          # MT/Ribo/HB percentages
    └── 3.Scatterplot.pdf       # Feature correlation

## 4. Key Parameters
Parameter	Default	Description
min.cells	5	Gene detected in ≥5 cells
min.features	300	Cell with ≥300 genes
percent_mito cutoff	<20%	Mitochondrial QC threshold
scale.factor	10,000	Normalization parameter

# 5. Dependencies
Package	Version	Purpose
Seurat	5.0	Single-cell analysis
tidyverse	2.0	Data wrangling
clustree	0.5	Cluster evaluation
cowplot	1.1	Plot arrangement

# 6. Citation
Please cite:

1.Deng G, Liang X, Pan Y et al. Single-Cell Transcriptomic Analysis of Different Liver Fibrosis Models: Elucidating Molecular Distinctions and Commonalities

2.Seurat:

Satija R, et al. Nat Biotechnol. 2015. doi:10.1038/nbt.3192

3.Original data papers:

GSE171904: [Yang W, He H, Wang T, Su N et al. Single-Cell Transcriptomic Analysis Reveals a Hepatic Stellate Cell-Activation Roadmap and Myofibroblast Origin During Liver Fibrosis in Mice. Hepatology 2021 Nov;74(5):2774-2790. PMID: 34089528]

GSE199638: [Mooring M, Yeung GA, Luukkonen P, Liu S et al. Hepatocyte CYR61 polarizes profibrotic macrophages to orchestrate NASH fibrosis. Sci Transl Med 2023 Sep 27;15(715):eade3157. PMID: 37756381]

GSE221481: [Zhang L, Zhao C, Dai W, Tong H et al. Disruption of cholangiocyte-B cell crosstalk by blocking the CXCL12-CXCR4 axis alleviates liver fibrosis. Cell Mol Life Sci 2023 Nov 27;80(12):379. PMID: 38010435]

# 7. FAQ
❓ How to modify filtering thresholds?
Edit these lines in the script:
r
sce.all <- subset(sce.all, 
                 subset = nFeature_RNA > 200 &
                          nFeature_RNA < 5000 &
                          percent_mito < 20)
❓ Where to add new datasets?
Duplicate one of the existing processing blocks (lines 10-85 for GSE171904) and adjust sample names.

8. License
MIT License © [Guofei Deng] 2024




