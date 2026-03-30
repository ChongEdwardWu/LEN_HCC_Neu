### Step 1.5: pySCENIC import/export
### Export merged RNA counts to LOOM and merge pySCENIC results back into Seurat

## ----------------------------- Section 0: Setup -----------------------------
rm(list = ls()); gc()
# Optional: load a personal R profile if you use one locally
if (file.exists("path_to_.radian_profile")) source("path_to_.radian_profile")
# .libPaths()

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(SingleCellExperiment)
  library(SCENIC)
  library(SCopeLoomR)
  library(data.table)
})

set.seed(123)

# Project paths
WORKDIR     <- "path_to_data"
RESULTS_DIR <- file.path(WORKDIR, "results")
FIG_DIR     <- file.path(WORKDIR, "figures")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
if (!dir.exists(FIG_DIR))     dir.create(FIG_DIR, recursive = TRUE)

# Input Seurat object (after integration/clustering)
SEURAT_RDS <- file.path(RESULTS_DIR, "01_Len_Integration.rds")
stopifnot(file.exists(SEURAT_RDS))

# pySCENIC outputs (produced by your Python script)
SCENIC_BASEDIR <- "path_to_pyscenic_results"
SCENIC_RES     <- file.path(SCENIC_BASEDIR, "results")
DATASET_ID     <- "01_pySCENIC_counts_merged"

# File names produced by run_SCENIC.py
SCENIC_LOOM   <- file.path(SCENIC_RES, paste0(DATASET_ID, ".scenic.loom"))
SCENIC_MOTIFS <- file.path(SCENIC_RES, paste0(DATASET_ID, ".motifs.csv"))
SCENIC_AUC    <- file.path(SCENIC_RES, paste0(DATASET_ID, ".auc.csv"))
SCENIC_BIN    <- file.path(SCENIC_RES, paste0(DATASET_ID, ".bin.csv"))
SCENIC_THR    <- file.path(SCENIC_RES, paste0(DATASET_ID, ".thresholds.csv"))

# Loom to export for pySCENIC step-1 (produced here)
EXPORT_LOOM   <- file.path(RESULTS_DIR, paste0(DATASET_ID, ".loom"))

# Gene filtering threshold for export (proportion of cells with >0 UMI)
GENE_MIN_PROP <- 0.005   # 0.5%; change to 0.01 if you prefer 1%


## ---------------------- Section 1: Build merged counts ----------------------
message("Loading Seurat object: ", SEURAT_RDS)
seu <- readRDS(SEURAT_RDS)

DefaultAssay(seu) <- "RNA"
# 1) Merge per-sample layers inside the RNA assay 
seu <- JoinLayers(seu)
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)

# 2) Extract the unified raw counts matrix (sparse; columns already aligned).
exprMat <- GetAssayData(seu, layer = "counts")

# 3) Ensure unique gene names (drop duplicated rows if present).
if (any(duplicated(rownames(exprMat)))) {
  dup_n <- sum(duplicated(rownames(exprMat)))
  message("Removing duplicated genes: ", dup_n)
  exprMat <- exprMat[!duplicated(rownames(exprMat)), , drop = FALSE]
}

# 4) Light gene filtering for pySCENIC (keep genes detected in >= GENE_MIN_PROP of cells).
stopifnot(exists("GENE_MIN_PROP"))
keep <- Matrix::rowSums(exprMat > 0) >= ceiling(GENE_MIN_PROP * ncol(exprMat))
exprMat_f <- exprMat[keep, , drop = FALSE]
message("Merged counts after filtering: ", paste(dim(exprMat_f), collapse = " x "))

# 5) Export to LOOM for pySCENIC.
if (file.exists(EXPORT_LOOM)) file.remove(EXPORT_LOOM)
loom <- build_loom(file.name = EXPORT_LOOM, dgem = exprMat_f)
close_loom(loom)
message("LOOM exported for pySCENIC: ", EXPORT_LOOM)


## ---------------- Section 2: Load pySCENIC results --------------------------
# You can read CSV or LOOM; both are supported here. We will then apply a unified
# "safe naming" mapping to regulon names below.

# --- Helper: standardize regulon names to *_pos / *_neg and ensure uniqueness ---
sanitize_regulon_names <- function(x) {
  x <- trimws(as.character(x))
  # Match only trailing (+) / (-)
  x <- sub("\\s*\\(\\+\\)\\s*$", "_pos", x, perl = TRUE)
  x <- sub("\\s*\\(-\\)\\s*$", "_neg", x, perl = TRUE)
  # Replace consecutive whitespace with one underscore to avoid odd separators
  x <- gsub("[[:space:]]+", "_", x, perl = TRUE)
  x
}

# 2.1 Read LOOM (optional) if present
regulons_incidMat <- regulons <- regulonAUC_loom <- NULL
if (file.exists(SCENIC_LOOM)) {
  message("Reading regulons & AUC from LOOM: ", SCENIC_LOOM)
  lm <- open_loom(SCENIC_LOOM)
  regulons_incidMat <- get_regulons(lm, column.attr.name = "Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC_loom <- get_regulons_AUC(lm, column.attr.name = "RegulonsAUC")
  close_loom(lm)
} else {
  message("LOOM not found (skip LOOM import): ", SCENIC_LOOM)
}

# 2.2 Read motif enrichment and AUC/Bin from CSV (robust and easy to merge)
files2check <- c(SCENIC_MOTIFS, SCENIC_AUC, SCENIC_BIN, SCENIC_THR)
if (!all(file.exists(files2check))) {
  stop("Missing SCENIC result files:\n", paste(files2check[!file.exists(files2check)], collapse = "\n"))
}

# Motif enrichment table (not strongly bound to regulon naming; no changes needed)
motifEnrichment <- fread(SCENIC_MOTIFS, header = TRUE, skip = 1)[-1, ]
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")

# AUC: pySCENIC CSV is "cell x regulon"; transpose to "regulon x cell"
regulonAUC2 <- fread(SCENIC_AUC, header = TRUE) |>
  column_to_rownames("Cell") |>
  t()

# Bin (binary activity): same handling
regulonBin <- fread(SCENIC_BIN, header = TRUE) |>
  column_to_rownames("Cell") |>
  t()

# ---- Apply unified "safe naming" and ensure AUC & Bin consistency ----
# 1) Build a mapping (orig -> safe) from original AUC regulon names
orig_auc_names  <- rownames(regulonAUC2)
safe_auc_names  <- sanitize_regulon_names(orig_auc_names)

# If duplicates remain (e.g., multiple variants of the same TF), make names unique and warn
if (any(duplicated(safe_auc_names))) {
  dups <- unique(safe_auc_names[duplicated(safe_auc_names)])
  warning(sprintf(
    "Detected %d duplicated regulon names after sanitization (e.g., %s). Appending suffix _vN to make them unique.",
    length(dups), paste(head(dups, 3), collapse = ", ")
  ))
  safe_auc_names <- make.unique(safe_auc_names, sep = "_v")
}
name_map <- setNames(safe_auc_names, orig_auc_names) # mapping: original -> safe

# 2) Apply to AUC
rownames(regulonAUC2) <- unname(name_map[rownames(regulonAUC2)])
regulonAUC2 <- Matrix(regulonAUC2, sparse = TRUE)

# 3) Apply to Bin (first sanitize similarly, then enforce alignment using AUC map)
orig_bin_names <- rownames(regulonBin)
tmp_bin_names  <- sanitize_regulon_names(orig_bin_names)
# Prefer the AUC-based mapping for consistency; if new names exist in Bin but not AUC, fall back
mapped_bin_names <- name_map[orig_bin_names]
need_fallback    <- is.na(mapped_bin_names)
if (any(need_fallback)) {
  warning(sprintf(
    "Found %d regulons in BIN not present in AUC; keeping their sanitized names.",
    sum(need_fallback)
  ))
  mapped_bin_names[need_fallback] <- tmp_bin_names[need_fallback]
}
rownames(regulonBin) <- unname(mapped_bin_names)
regulonBin <- Matrix(regulonBin, sparse = TRUE)

# 4) (Optional) Also rename objects loaded from LOOM to keep naming consistent
if (!is.null(regulons)) {
  # 'regulons' is a list with regulon names as element names
  n_old <- names(regulons)
  n_new <- ifelse(n_old %in% names(name_map), unname(name_map[n_old]), sanitize_regulon_names(n_old))
  # If duplicates remain, make unique (use the same rule as AUC)
  if (any(duplicated(n_new))) n_new <- make.unique(n_new, sep = "_v")
  names(regulons) <- n_new
}
if (!is.null(regulons_incidMat)) {
  # 'regulons_incidMat' is usually a sparse matrix (regulon x gene) or similar
  rn_old <- rownames(regulons_incidMat)
  rn_new <- ifelse(rn_old %in% names(name_map), unname(name_map[rn_old]), sanitize_regulon_names(rn_old))
  if (any(duplicated(rn_new))) rn_new <- make.unique(rn_new, sep = "_v")
  rownames(regulons_incidMat) <- rn_new
}
if (!is.null(regulonAUC_loom)) {
  # If 'regulonAUC_loom' is a matrix, rename its rownames
  rn_old <- rownames(regulonAUC_loom)
  rn_new <- ifelse(rn_old %in% names(name_map), unname(name_map[rn_old]), sanitize_regulon_names(rn_old))
  if (any(duplicated(rn_new))) rn_new <- make.unique(rn_new, sep = "_v")
  rownames(regulonAUC_loom) <- rn_new
}

# 5) Thresholds table (optional): keep original column and add 'SafeName' for cross-reference
regulonThr <- fread(SCENIC_THR, header = TRUE)
if (!"Regulon" %in% colnames(regulonThr)) {
  # Column names can vary across versions; try to guess a proper column
  cand <- intersect(colnames(regulonThr), c("regulon", "Name", "name", "TF"))
  if (length(cand)) {
    setnames(regulonThr, cand[1], "Regulon")
  } else {
    warning("Cannot find 'Regulon' column in thresholds table; skipping SafeName addition.")
  }
}
if ("Regulon" %in% colnames(regulonThr)) {
  regulonThr$SafeName <- ifelse(
    regulonThr$Regulon %in% names(name_map),
    unname(name_map[regulonThr$Regulon]),
    sanitize_regulon_names(regulonThr$Regulon)
  )
}

## ---------------- Section 3: Merge results into Seurat ----------------------
# Reload Seurat to ensure a clean object
seu <- readRDS(SEURAT_RDS)

# Ensure identical cell sets and reorder columns to match Seurat
stopifnot(all(colnames(seu) %in% colnames(regulonAUC2)))
regulonAUC2 <- regulonAUC2[, colnames(seu), drop = FALSE]

stopifnot(all(colnames(seu) %in% colnames(regulonBin)))
regulonBin  <- regulonBin[,  colnames(seu), drop = FALSE]

# Final robustness checks: unique rownames and AUC/Bin name consistency
stopifnot(!any(duplicated(rownames(regulonAUC2))))
stopifnot(!any(duplicated(rownames(regulonBin))))
stopifnot(identical(sort(rownames(regulonAUC2)), sort(rownames(regulonBin))))

# Save an "original -> safe" mapping table for traceability
regulon_name_map <- tibble(
  regulon_orig = names(name_map),
  regulon_safe = unname(name_map)
)

# Add to Seurat
seu[["AUC"]] <- CreateAssayObject(data = regulonAUC2)
seu[["Bin"]] <- CreateAssayObject(data = regulonBin)
DefaultAssay(seu) <- "SCT"

# Save related objects and the name mapping together
save(regulons_incidMat, regulons, regulonAUC_loom,
     motifEnrichment, regulonAUC2, regulonBin, regulonThr, regulon_name_map,
     file = file.path(RESULTS_DIR, "01.5_LEN_pySCENIC_results.RData"))

# Save the updated Seurat object
SAVE_RDS <- file.path(RESULTS_DIR, "01.5_LEN_pySCENIC_results.rds")
saveRDS(seu, SAVE_RDS)
message("Updated Seurat saved: ", SAVE_RDS)

# Session info
writeLines(capture.output(sessionInfo()),
           file.path(WORKDIR, "session_info_01.5_LEN_pySCENIC_results.txt"))
