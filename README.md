# LEN_HCC_Neu

Analysis code for single-cell RNA sequencing, regulon, neutrophil-state, and lineage-resolved treatment-response analyses in a mouse hepatocellular carcinoma model treated with lenvatinib.

The repository contains path-sanitized code for preprocessing, quality control, Seurat integration, immune-cell annotation, neutrophil reference mapping, pySCENIC analysis, differential expression, enrichment analysis, and figure generation. Raw data and large external reference objects are not included.

## Repository structure

```text
LEN_HCC_Neu/
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
│       ├── Step5_Len_LENvsCTR_Figures.R
│       ├── Step6_Len_Neu_query_native_display.R
│       └── Step7_Len_CD45_lineage_remodeling.R
├── python/
│   └── run_SCENIC_mouse.py
├── shell/
│   └── run_cellbender_batch.sh
├── env/
│   └── .gitkeep
└── README.md
```

## Analysis workflow

### Core workflow

1. `shell/run_cellbender_batch.sh`
   - Runs `cellbender remove-background` on Cell Ranger raw-feature matrices.

2. `R/QC/scRNA_QC_loop_CB_step1.R`
   - Reads CellBender-filtered matrices, calculates sample-level quality-control metrics, performs exploratory clustering, and calls doublets with `scDblFinder`.

3. `R/QC/scRNA_QC_loop_CB_step2.R`
   - Applies final filtering, harmonizes gene symbols, calculates cell-cycle scores, and creates sample-level Seurat objects with `SCTransform`.

4. `R/Analysis/Step1_Len_Data_Integration.R`
   - Performs multi-sample SCT-RPCA integration, dimensionality reduction, and graph-based clustering.

5. `R/Analysis/Step1_5_Len_pySCENIC.R` and `python/run_SCENIC_mouse.py`
   - Export merged RNA counts, run pySCENIC, and add regulon activity to the Seurat object.

6. `R/Analysis/Step2_Len_Clustering.R`
   - Assigns broad immune-cell annotations and summarizes cluster composition.

7. `R/Analysis/Step3_Len_Neu_sub_clustering.R`
   - Subsets neutrophils, maps discrete Ng et al. maturation states, and prepares the neutrophil object used by downstream analyses.

8. `R/Analysis/Step4_Len_Neu_DEG.R`
   - Performs neutrophil differential-expression, regulon-activity, score, and gene-set enrichment analyses.

9. `R/Analysis/Step5_Len_LENvsCTR_Figures.R`
   - Generates neutrophil-focused figure panels.

10. `R/Analysis/Step6_Len_Neu_query_native_display.R`
    - Recalculates a query-native neutrophil UMAP, performs de novo five-cluster analysis, retains upstream Ng state labels, transfers tumour-derived Xue et al. mouse neutrophil states, and calculates the Xue core-state maturation score.

11. `R/Analysis/Step7_Len_CD45_lineage_remodeling.R`
    - Performs matched-background differential-expression analyses across eight immune lineages, identifies genes with concordant attenuation under PERK inhibition and myeloid Eif2ak3 deletion, runs GO Biological Process enrichment, and generates the myeloid and lymphoid Venn/GO displays.

## Main software dependencies

### R

- `Seurat`, `SeuratObject`, `SeuratWrappers`
- `SingleCellExperiment`, `scater`, `scran`, `scuttle`, `DropletUtils`
- `scDblFinder`, `bluster`, `BiocSingular`, `BiocParallel`
- `SingleR`, `UCell`, `destiny`, `presto`
- `SCENIC`, `SCopeLoomR`
- `msigdbr`, `fgsea`, `clusterProfiler`, `org.Mm.eg.db`, `AnnotationDbi`, `GOSemSim`
- `tidyverse`, `Matrix`, `openxlsx2`, `patchwork`, `ggridges`, `ggVennDiagram`, `ragg`

### Python and command-line tools

- `cellbender`
- `pyscenic`
- `loompy`
- `pandas`

Step 6 requires Seurat 5.1 or later.

## External inputs

- Cell Ranger and CellBender outputs
- pySCENIC motif annotations, transcription-factor lists, and ranking databases
- a processed Ng et al. neutrophil reference object
- the Ng et al. supplementary gene list used by the original Step 3 workflow
- a processed Xue et al. mouse liver-tumour reference object
- the annotated CD45-positive Seurat object required by Step 7
- any additional gene-set resources referenced by the core workflow

Replace placeholder paths in the core scripts and supply Step 6 and Step 7 paths through their command-line arguments.

## Reproducibility and inference boundaries

- Random seeds are fixed in the analysis scripts.
- Step 6 verifies that reference mapping and de novo clustering do not alter the independently calculated query-native UMAP coordinates.
- Step 7 records input checksums, parameters, session information, cell counts, tested-gene universes, and fold-change quality-control results.
- Each experimental condition is represented by one pooled 10x library. Cell-level tests and adjusted P values are exploratory distributional summaries and do not provide treatment-level biological-replicate inference.
- Cell-state proportions and score distributions are descriptive. The Step 7 attenuation criterion is not a formal interaction or rescue test.

## Repository scope

This repository includes the computational analysis code used for the single-cell and neutrophil-state analyses. It does not include raw sequencing data, large external references, manuscript source files, or wet-lab analysis code unrelated to the single-cell workflow.

Refer to the manuscript Data and Code Availability section for accession numbers and data-sharing details.

