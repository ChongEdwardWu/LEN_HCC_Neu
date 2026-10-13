# LEN_HCC_Neu

Analysis code for mouse single-cell and bulk RNA sequencing and human nCounter immune-signature analyses in hepatocellular carcinoma treated with lenvatinib.

The repository contains path-sanitized code for preprocessing, quality control, Seurat integration, immune-cell annotation, neutrophil reference mapping, pySCENIC analysis, differential expression, enrichment analysis, and figure generation. Raw data, patient-level inputs, and large external reference objects are not included.

## Repository structure

```text
LEN_HCC_Neu/
├── R/
│   ├── QC/
│   │   ├── scRNA_QC_loop_CB_step1.R
│   │   └── scRNA_QC_loop_CB_step2.R
│   ├── BulkRNA/
│   │   └── Step01_Neu_RNA_DEG.R
│   ├── Human/
│   │   ├── Step1_Human_nCounter_signature_scores.R
│   │   ├── Step2_Human_PD_delta_heatmaps.R
│   │   └── Step3_Human_PD_signature_ranking.R
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

### Human nCounter workflow

The three scripts in `R/Human/` cover Figures S1A, 1A, and S1B, respectively. Run them from the repository root with `readxl`, `ggplot2`, and `pheatmap` installed:

```bash
Rscript R/Human/Step1_Human_nCounter_signature_scores.R "data/human/nCounter_Normalization-1.xlsx" "data/human/sample_metadata.csv" "results/human/signatures"
Rscript R/Human/Step2_Human_PD_delta_heatmaps.R "data/human/patient annotation-1.xlsx" "results/human/signatures/signature_scores_and_z_by_sample.csv" "results/human/paired_changes"
Rscript R/Human/Step3_Human_PD_signature_ranking.R "data/human/S1B_STable_PD_delta_stats_ALLsignatures.csv" "results/human/ranking"
```

- **Inputs:** normalized nCounter workbook (`NC_9785_N1_Normalization` sheet; sample IDs in row 3, samples from column 8, genes from row 5); sample metadata with `sample_id`, `patient_id`, and `timepoint` (`baseline`, `week4`, or `PD`); clinical workbook (`Patient annotation` sheet; `FL_ID` and `mRECIST`). A final argument can override each sheet name. The supplied sample-score table can also serve as metadata: only its three metadata columns are read. Metadata must match all expression samples, and RCC-header/FL-ID column selections must agree.
- **Scores and paired changes:** each signature is the mean of marker-level `log2(normalized expression + 1)`, not a cell fraction. Z scores use the mean and sample SD of all baseline samples. Step 2 calculates within-patient `PD biopsy Z - baseline Z`, excludes week-4 biopsies, and summarizes by mRECIST response (`PD` versus `CR/PR/SD`, grouped as `non-PD`). The biopsy timepoint and response group are separate variables. `CXCR2/CSF3R` is consistently named `Neutrophil`.
- **Marker plot:** color is `log2(mean normalized expression + 1)`, distinct from the signature-score calculation. Dot size uses the workbook percentage above threshold; only if that column is absent does it use the fraction of expression values greater than zero. Missing or invalid threshold percentages are not silently replaced.
- **Ranking:** Step 3 reads the supplied summary statistics, including `signature`, `diff_median`, `prop_up_PD`, `prop_up_nonPD`, `r_rank_biserial`, `ci_med_low`, and `ci_med_high`. It retains the descriptive formula `0.5*z(diff_median) + 0.2*z(prop_up_PD - prop_up_nonPD) + 0.3*z(sign(diff_median)*abs(r_rank_biserial))`, with standardization across signatures. The x-axis shows the median difference; signatures are ordered by the composite score. Confidence intervals and inferential statistics are read from the input table, not recalculated by these scripts.

Each script writes figures and supporting CSV tables to a new or empty output directory and records `sessionInfo.txt`. No raw nCounter normalization or survival analysis is performed. Keep patient-level inputs and outputs outside version control.

## Main software dependencies

### R

- `Seurat`, `SeuratObject`, `SeuratWrappers`
- `SingleCellExperiment`, `scater`, `scran`, `scuttle`, `DropletUtils`
- `scDblFinder`, `bluster`, `BiocSingular`, `BiocParallel`
- `SingleR`, `UCell`, `destiny`, `presto`
- `SCENIC`, `SCopeLoomR`
- `msigdbr`, `fgsea`, `clusterProfiler`, `org.Mm.eg.db`, `AnnotationDbi`, `GOSemSim`
- `tidyverse`, `Matrix`, `openxlsx2`, `patchwork`, `ggridges`, `ggVennDiagram`, `ragg`
- `readxl`, `ggplot2`, `pheatmap` (human nCounter workflow)

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
- In the mouse scRNA-seq dataset, each experimental condition is represented by one pooled 10x library. Cell-level tests and adjusted P values are exploratory distributional summaries and do not provide treatment-level biological-replicate inference.
- Cell-state proportions and score distributions are descriptive. The Step 7 attenuation criterion is not a formal interaction or rescue test.

## Repository scope

This repository includes mouse single-cell, neutrophil-state, and bulk RNA-seq code, together with human nCounter immune-signature scoring and visualization. It does not include raw sequencing data, patient-level data, large external references, or manuscript source files.

Refer to the manuscript Data and Code Availability section for accession numbers and data-sharing details.

