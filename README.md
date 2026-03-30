# LEN_HCC_Neu

Public code repository for the Lenvatinib/PERK single-cell RNA-seq analysis workflow, with a focus on neutrophil-state analysis in HCC.

## Repository layout

- `R/Analysis/`: main R analysis scripts from integration to DEG analysis and figure generation.
- `R/QC/`: sample-level preprocessing and QC scripts.
- `python/`: helper Python scripts used for pySCENIC processing.
- `shell/`: helper shell scripts for batch preprocessing steps.
- `env/`: placeholder directory for environment files.

## Included R workflow

- `Step1_Len_Data_Integration.R`
- `Step1_5_Len_pySCENIC.R`
- `Step2_Len_Clustering.R`
- `Step3_Len_Neu_sub_clustering.R`
- `Step4_Len_Neu_DEG.R`
- `Step5_Len_LENvsCTR_Figures.R`
- `R/QC/scRNA_QC_loop_CB_step1.R`
- `R/QC/scRNA_QC_loop_CB_step2.R`

These are sanitized public versions:

- absolute personal paths were replaced with placeholders
- Chinese comments were translated into English
- date suffixes were removed from output file names
- section numbering was normalized
- personal logging, email notifications, and clearly non-essential exploratory code were removed where appropriate

## Placeholder paths to update

Before running the scripts, replace placeholders such as:

- `path_to_project_root`
- `path_to_data`
- `path_to_scenic_reference`
- `path_to_bulk_deg_xlsx`
- `path_to_ng_reference_rds`
- `path_to_ng_table_s1_xlsx`
- `path_to_figs2a_moesm4_xlsx`
- `path_to_immgen_reference_rdata`
- `path_to_azimuth_human_reference_dir`
- `path_to_azimuth_mouse_reference_dir`
- `path_to_cellranger_outputs`
- `path_to_cellbender_output`
- `path_to_conda_sh`

## Notes

- `Step3_Len_Neu_sub_clustering.R` expects the public Step 2 output `results/02_Annotation_LSK.rds`.
- `Step4_Len_Neu_DEG.R` writes Excel outputs into `results/Step4_GroupDE_GSEA/`.
- `Step5_Len_LENvsCTR_Figures.R` writes final figures into `figures/final/`.
- `python/run_SCENIC_mouse.py` and `shell/run_cellbender_batch.sh` are included as helper scripts and may still need project-specific tuning before execution.
