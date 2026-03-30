# LEN_HCC_Neu

Code repository for the manuscript:

**A therapy-induced PERK–ATF4 checkpoint in neutrophils drives immune exclusion and resistance to lenvatinib in hepatocellular carcinoma**

Repository URL: <https://github.com/ChongEdwardWu/LEN_HCC_Neu>

## Overview

This repository contains the public analysis scripts used for the single-cell RNA-seq, pySCENIC, neutrophil-state mapping, differential-expression, and figure-generation workflows in the manuscript. The code is organized to follow the order in which the analyses were performed, from preprocessing and quality control to integration, neutrophil-focused downstream analysis, and final figure generation.

The repository is intended as a **scripts-focused companion** to the manuscript. Large raw data files, intermediate objects, and external reference resources are **not** stored here.

## Repository structure

```text
LEN_HCC_Neu/
├── README.md
├── .gitignore
├── R/
│   ├── QC/
│   │   ├── scRNA_QC_loop_CB_step1.R
│   │   └── scRNA_QC_loop_CB_step2.R
│   └── Analysis/
│       ├── Step1_Len_Data_Integration.R
│       ├── Step1_5_Len_pySCENIC.R
│       ├── Step2_Len_Clustering.R
│       ├── Step3_Len_Neu_sub_clustering.R
│       ├── Step4_Len_Neu_DEG.R
│       └── Step5_Len_LENvsCTR_Figures.R
├── python/
│   └── run_SCENIC_mouse.py
├── shell/
│   └── run_cellbender_batch.sh
└── env/
    └── .gitkeep
```

## Workflow summary

### 1. Preprocessing and quality control

- **`shell/run_cellbender_batch.sh`**  
  Runs `cellbender remove-background` on Cell Ranger `raw_feature_bc_matrix` directories.

- **`R/QC/scRNA_QC_loop_CB_step1.R`**  
  Imports CellBender-corrected `.h5` files, performs initial normalization, feature selection, clustering, doublet detection, and preliminary QC visualization.

- **`R/QC/scRNA_QC_loop_CB_step2.R`**  
  Applies final filtering, creates Seurat objects, performs cell-cycle scoring and SCTransform normalization, and prepares per-sample objects for integration.

### 2. Integration and global clustering

- **`R/Analysis/Step1_Len_Data_Integration.R`**  
  Integrates QC-passed samples using Seurat SCT-RPCA, performs PCA, UMAP, and graph-based clustering, and saves the integrated object.

### 3. Regulon analysis with pySCENIC

- **`R/Analysis/Step1_5_Len_pySCENIC.R`**  
  Exports merged RNA counts to loom format for pySCENIC and imports pySCENIC outputs back into Seurat as AUC and binary regulon assays.

- **`python/run_SCENIC_mouse.py`**  
  Runs pySCENIC on the exported loom file using:
  - GRNBoost2 for gene regulatory network inference
  - cisTarget motif enrichment
  - AUCell regulon scoring
  - binarization of regulon activity

### 4. Neutrophil-focused downstream analysis

- **`R/Analysis/Step2_Len_Clustering.R`**  
  Performs cluster-level annotation of the integrated CD45+ immune-cell dataset.

- **`R/Analysis/Step3_Len_Neu_sub_clustering.R`**  
  Subsets neutrophils, maps them to the Ng et al. neutrophil reference atlas, computes neutrophil maturation scores, and prepares neutrophil state objects for downstream analyses.

- **`R/Analysis/Step4_Len_Neu_DEG.R`**  
  Performs differential expression, regulon comparisons, score-based comparisons, and GSEA across treatment groups within neutrophils.

### 5. Figure generation

- **`R/Analysis/Step5_Len_LENvsCTR_Figures.R`**  
  Generates figure panels and summary plots used in the manuscript.

## Suggested execution order

The scripts are not fully wrapped into an automated pipeline. A typical execution order is:

1. `shell/run_cellbender_batch.sh`
2. `R/QC/scRNA_QC_loop_CB_step1.R`
3. `R/QC/scRNA_QC_loop_CB_step2.R`
4. `R/Analysis/Step1_Len_Data_Integration.R`
5. `R/Analysis/Step1_5_Len_pySCENIC.R` (export loom)
6. `python/run_SCENIC_mouse.py`
7. `R/Analysis/Step1_5_Len_pySCENIC.R` (import pySCENIC results)
8. `R/Analysis/Step2_Len_Clustering.R`
9. `R/Analysis/Step3_Len_Neu_sub_clustering.R`
10. `R/Analysis/Step4_Len_Neu_DEG.R`
11. `R/Analysis/Step5_Len_LENvsCTR_Figures.R`

## External inputs required

This public repository does **not** include large or third-party input files. To run the workflow, you will need to provide the following resources locally:

1. **Cell Ranger outputs** for each sample
2. **CellBender-corrected `.h5` files**
3. **pySCENIC reference resources**, including:
   - motif ranking databases
   - motif annotation tables
   - transcription factor lists
4. **Reference objects and supplementary resources** used by downstream scripts, including:
   - Ng et al. neutrophil reference Seurat object
   - Ng et al. Table S1 gene list for maturation scoring
   - figure-related gene-set spreadsheets used in the manuscript
   - optional ImmGen / Azimuth reference files if you enable those sections

## Placeholder paths to update

The public scripts are sanitized and still contain placeholders that must be edited before use. Common examples include:

- `path_to_project_root`
- `path_to_data`
- `path_to_cellranger_outputs`
- `path_to_cellbender_output`
- `path_to_conda_sh`
- `path_to_scenic_reference`
- `path_to_pyscenic_results`
- `path_to_ng_reference_rds`
- `path_to_ng_table_s1_xlsx`
- `path_to_figs2a_moesm4_xlsx`
- `path_to_bulk_deg_xlsx`
- `path_to_immgen_reference_rdata`
- `path_to_azimuth_human_reference_dir`
- `path_to_azimuth_mouse_reference_dir`

## Software environment

The analysis was developed in a Linux environment using:

- **R 4.3.3**
- **Seurat 5.1.0**
- **scDblFinder 1.14.0**
- **CellBender 0.3.0**
- **pySCENIC 0.12.1**

Additional R package dependencies are declared within the individual scripts.

## Data availability

According to the manuscript, the sequencing datasets are deposited in the Genome Sequence Archive (GSA):

- **scRNA-seq**: CRA038908
- **bulk RNA-seq**: CRA038902

Please refer to the manuscript for current release status and accession details.

## Reproducibility notes

- These are **sanitized public versions** of the original scripts.
- Absolute local paths and non-essential personal logging have been removed.
- Some scripts expect intermediate objects generated by earlier steps.
- Large intermediate files such as `.rds`, `.loom`, `.h5`, and figure caches are intentionally not versioned in this repository.
- Random seeds are fixed in key steps where explicitly coded.

## Contact

For questions about the code or the manuscript, please contact the corresponding authors listed in the manuscript.

