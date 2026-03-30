## =======================================================================
## Step4: Neutrophil group DE + GSEA + Reversal highlights
##   - RNA DE (FindMarkers) + GSEA (fgsea)
##   - Regulons: AUC / Bin DE
##   - Scores: UCell + NeutrophilMaturation_UCell (cell-level)
##   - AllNeu only (NO scopes)
##   - Output: per-comparison Excel + one Reversal highlight Excel
## =======================================================================

## -------------------------- Section 0: Preparation -----------------------
rm(list = ls())
gc()

# Optional: load a personal R profile if you use one locally
if (file.exists("path_to_.radian_profile")) source("path_to_.radian_profile")

workdir <- "path_to_data"
setwd(workdir)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(openxlsx2)
  library(future)
  library(msigdbr)
  library(fgsea)
})

set.seed(123)

# parallel (Linux OK). If unstable, switch to sequential.
nworkers <- 8
future::plan("multicore", workers = nworkers)
# future::plan("sequential")

# I/O
# Make sure this input RDS contains the UCell and FigS2A score columns.
in_rds <- file.path(workdir, "results/03_Len_Neu_NgMapping.rds")
stopifnot(file.exists(in_rds))

out_dir <- file.path(workdir, "results", "Step4_GroupDE_GSEA")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


## -------------------------- Section 1: Load + basic checks ----------------
seu <- readRDS(in_rds)

stopifnot("group" %in% colnames(seu@meta.data))
stopifnot("Ng_stage" %in% colnames(seu@meta.data))

group_levels <- c("CTR","CTRnAMG","CTRnKO","LEN","LENnAMG","LENnKO")
seu$group <- factor(as.character(seu$group), levels = group_levels)

message("Cells: ", ncol(seu), " | Genes: ", nrow(seu))
print(table(seu$group, useNA = "ifany"))
print(table(seu$Ng_stage, useNA = "ifany"))

# keep only neutrophils if column exists
if ("CellType_l1" %in% colnames(seu@meta.data)) {
  seu <- subset(seu, subset = CellType_l1 == "Neu")
  message("After subset Neu: cells=", ncol(seu))
}

stopifnot("RNA" %in% names(seu@assays))
stopifnot("AUC" %in% names(seu@assays))
stopifnot("Bin" %in% names(seu@assays))


## -------------------------- Section 2: Feature blacklist (RNA) ------------
build_bad_features <- function(gene_names) {
  hist_genes    <- grep("^Hist", gene_names, ignore.case = TRUE, value = TRUE)
  hb_genes      <- grep("^Hb[ab]-|^HB[^(P)]", gene_names, value = TRUE)
  mt_genes      <- grep("^mt-|^MT-|MTRNR2L|Mtrnr2l|^Mtmr", gene_names, ignore.case = TRUE, value = TRUE)
  rps_genes     <- grep("^Rp[sl]|^RP[SL]", gene_names, ignore.case = TRUE, value = TRUE)
  rik_genes     <- grep("^.*Rik$|^Rik", gene_names, ignore.case = TRUE, value = TRUE)
  pseudo_genes  <- grep("-rs|-ps", gene_names, ignore.case = TRUE, value = TRUE)
  mir_genes     <- grep("^Mir", gene_names, ignore.case = TRUE, value = TRUE)
  gencode_genes <- grep("^Gm|ENSMUSG", gene_names, ignore.case = TRUE, value = TRUE)

  unique(c(hist_genes, hb_genes, mt_genes, rps_genes, rik_genes, pseudo_genes, mir_genes, gencode_genes))
}

bad_features <- build_bad_features(rownames(seu))
features_rna <- setdiff(rownames(seu), bad_features)

message("Bad features: ", length(bad_features))
message("RNA features used: ", length(features_rna))


## -------------------------- Section 3: MSigDB pathways for GSEA -----------
# Pre-load once (avoid repeated msigdbr calls)
build_pathways <- function(category, subcategory = NULL, species = "Mus musculus") {
  df <- if (is.null(subcategory)) {
    msigdbr::msigdbr(species = species, category = category)
  } else {
    msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)
  }
  # gene_symbol should exist
  pathways <- split(df$gene_symbol, df$gs_name)
  pathways <- lapply(pathways, unique)
  pathways
}

pathways_REACTOME <- build_pathways("C2", "CP:REACTOME")
pathways_KEGG     <- build_pathways("C2", "CP:KEGG")
pathways_GOBP     <- build_pathways("C5", "GO:BP")
pathways_HALLMARK <- build_pathways("H")

message("MSigDB loaded: Reactome=", length(pathways_REACTOME),
        " | KEGG=", length(pathways_KEGG),
        " | GO:BP=", length(pathways_GOBP),
        " | Hallmark=", length(pathways_HALLMARK))


## -------------------------- Section 4: Helpers (Excel/DE/GSEA/Scores) -----
safe_sheet <- function(x, max_len = 31) {
  x <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", x)
  x <- gsub("\\s+", "_", x)
  substr(x, 1, max_len)
}

write_excel_openxlsx2 <- function(sheets, out_xlsx) {
  # drop NULL
  keep <- !vapply(sheets, is.null, logical(1))
  sheets <- sheets[keep]

  # also drop 0-row data.frames (optional)
  keep2 <- vapply(sheets, function(x) {
    if (is.data.frame(x)) nrow(x) > 0 else TRUE
  }, logical(1))
  sheets <- sheets[keep2]

  nm <- names(sheets)
  nm2 <- make.unique(vapply(nm, safe_sheet, character(1)))
  names(sheets) <- nm2

  if ("write_xlsx" %in% getNamespaceExports("openxlsx2")) {
    openxlsx2::write_xlsx(sheets, file = out_xlsx)
  } else {
    wb <- openxlsx2::wb_workbook()
    for (sh in names(sheets)) {
      wb$add_worksheet(sh)
      wb$add_data(sh, x = as.data.frame(sheets[[sh]]))
    }
    wb$save(out_xlsx, overwrite = TRUE)
  }
  message("Excel saved -> ", out_xlsx)
}

detect_fc_col <- function(df) {
  cand <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  hit <- cand[cand %in% colnames(df)][1]
  if (is.na(hit)) stop("Cannot find FC column in df. Columns: ", paste(colnames(df), collapse = ", "))
  hit
}

run_findmarkers <- function(seu_obj,
                            assay,
                            ident1, ident2,
                            features = NULL,
                            min_pct = 0.1,
                            logfc = 0.25,
                            test_use = "wilcox",
                            min_cells = 50,
                            return_thresh = 1) {

  tab <- table(Seurat::Idents(seu_obj))
  if (!all(c(ident1, ident2) %in% names(tab))) return(NULL)

  if (tab[[ident1]] < min_cells || tab[[ident2]] < min_cells) {
    message("Skip FindMarkers: ", ident1, "=", tab[[ident1]], " ; ", ident2, "=", tab[[ident2]])
    return(NULL)
  }

  args <- list(
    object         = seu_obj,
    assay          = assay,
    ident.1        = ident1,
    ident.2        = ident2,
    features       = features,
    min.pct        = min_pct,
    logfc.threshold= logfc,
    test.use       = test_use,
    only.pos       = FALSE,
    return.thresh  = return_thresh
  )

  # Some Seurat versions support densify and some do not.
  fm_formals <- names(formals(Seurat::FindMarkers))
  if ("densify" %in% fm_formals) {
    args$densify <- TRUE
  }

  out <- do.call(Seurat::FindMarkers, args)

  out <- out %>%
    tibble::rownames_to_column(var = "feature") %>%
    dplyr::arrange(p_val_adj, p_val)

  out
}

run_fgsea_from_rank <- function(rank_df,
                                pathways,
                                minSize = 10,
                                maxSize = 500,
                                nperm = 10000,
                                nproc = 1) {

  if (is.null(rank_df) || nrow(rank_df) == 0) return(NULL)

  fc_col <- detect_fc_col(rank_df)

  stats <- rank_df[[fc_col]]
  names(stats) <- rank_df$feature

  stats <- stats[!is.na(stats)]
  stats <- stats[!is.na(names(stats))]
  stats <- stats[!duplicated(names(stats))]

  stats <- sort(stats, decreasing = TRUE)
  if (length(stats) < 100) return(NULL)

  fg <- fgsea::fgsea(
    pathways = pathways,
    stats    = stats,
    minSize  = minSize,
    maxSize  = maxSize,
    nperm    = nperm,
    nproc    = nproc
  )

  fg <- fg %>%
    as.data.frame() %>%
    dplyr::as_tibble() %>%
    dplyr::arrange(padj, pval)

  # collapse leadingEdge
  if ("leadingEdge" %in% colnames(fg)) {
    fg$leadingEdge <- vapply(fg$leadingEdge, function(x) paste(x, collapse = ", "), character(1))
  }

  fg
}

score_compare_celllevel <- function(seu_obj,
                                   score_cols,
                                   ident1, ident2,
                                   group_col = "group",
                                   min_cells = 50) {

  md <- seu_obj@meta.data
  keep <- intersect(score_cols, colnames(md))
  if (!length(keep)) return(NULL)

  df <- md[, c(group_col, keep), drop = FALSE]
  df[[group_col]] <- as.character(df[[group_col]])

  n1 <- sum(df[[group_col]] == ident1)
  n2 <- sum(df[[group_col]] == ident2)
  if (n1 < min_cells || n2 < min_cells) {
    message("Skip score_compare (too few cells): ", ident1, "=", n1, " ; ", ident2, "=", n2)
    return(NULL)
  }

  res <- lapply(keep, function(sc) {
    x1 <- df[df[[group_col]] == ident1, sc]
    x2 <- df[df[[group_col]] == ident2, sc]
    if (length(x1) < 10 || length(x2) < 10) return(NULL)

    wt <- suppressWarnings(stats::wilcox.test(x1, x2, exact = FALSE))
    tibble::tibble(
      score = sc,
      n1 = length(x1), n2 = length(x2),
      mean1 = mean(x1, na.rm = TRUE),
      mean2 = mean(x2, na.rm = TRUE),
      median1 = stats::median(x1, na.rm = TRUE),
      median2 = stats::median(x2, na.rm = TRUE),
      delta_mean = mean1 - mean2,
      delta_median = median1 - median2,
      p = wt$p.value
    )
  }) %>% dplyr::bind_rows()

  if (nrow(res) == 0) return(NULL)

  res <- res %>%
    dplyr::mutate(p_adj = p.adjust(p, method = "BH")) %>%
    dplyr::arrange(p_adj, p)

  res
}


## -------------------------- Section 5: Define comparisons (treat vs ctrl) -
# Direction: ident.1 is the treatment group and ident.2 is the reference group.
comparisons <- list(
  LEN_vs_CTR      = c("LEN", "CTR"),
  LENnAMG_vs_LEN  = c("LENnAMG", "LEN"),
  LENnKO_vs_LEN   = c("LENnKO", "LEN"),
  CTRnAMG_vs_CTR  = c("CTRnAMG", "CTR"),
  CTRnKO_vs_CTR   = c("CTRnKO", "CTR")   # corrected comparison key
)

# AllNeu only
scope_name <- "AllNeu"

# Score columns: maturation + all *_UCell
score_cols <- c(
  "NeutrophilMaturation_UCell",
  grep("_UCell$", colnames(seu@meta.data), value = TRUE)
) %>% unique()

message("Score columns detected: ", length(score_cols))


## -------------------------- Section 6: Run DE + GSEA per comparison -------
min_cells_per_group <- 50

# store results in memory for reversal summary
res_all <- list()

for (cmp_nm in names(comparisons)) {

  g1 <- comparisons[[cmp_nm]][1]  # treatment
  g2 <- comparisons[[cmp_nm]][2]  # control/reference

  # check groups exist in this object
  if (!all(c(g1, g2) %in% unique(as.character(seu$group)))) {
    message("Skip comparison (group missing in object): ", cmp_nm, " [", g1, " vs ", g2, "]")
    next
  }

  seu_sub <- subset(seu, subset = group %in% c(g1, g2))
  seu_sub$group <- factor(as.character(seu_sub$group), levels = c(g1, g2))
  Seurat::Idents(seu_sub) <- "group"

  n1 <- sum(seu_sub$group == g1)
  n2 <- sum(seu_sub$group == g2)
  if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
    message("Skip comparison (too few cells): ", cmp_nm, " ", g1, "=", n1, " ; ", g2, "=", n2)
    next
  }

  message("Running: ", cmp_nm, " | cells=", ncol(seu_sub),
          " (", g1, "=", n1, " ; ", g2, "=", n2, ")")

  # QC tables
  qc_cells_group <- as.data.frame(table(seu_sub$group))
  colnames(qc_cells_group) <- c("group", "n_cells")

  qc_cells_stage <- as.data.frame(table(seu_sub$group, seu_sub$Ng_stage))
  colnames(qc_cells_stage) <- c("group", "Ng_stage", "n_cells")

  qc_cells_sample <- NULL
  if ("Sample" %in% colnames(seu_sub@meta.data)) {
    qc_cells_sample <- seu_sub@meta.data %>%
      dplyr::group_by(group, Sample) %>%
      dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(group, dplyr::desc(n_cells))
  }

  ## ---- RNA DE (strict DE table)
  Seurat::DefaultAssay(seu_sub) <- "RNA"
  deg_rna <- run_findmarkers(
    seu_obj       = seu_sub,
    assay         = "RNA",
    ident1        = g1,
    ident2        = g2,
    features      = features_rna,
    min_pct       = 0.1,
    logfc         = 0.25,
    test_use      = "wilcox",
    min_cells     = min_cells_per_group,
    return_thresh = 1
  )

  ## ---- RNA rank table for GSEA (no logFC filter; include all genes)
  deg_rna_rank <- run_findmarkers(
    seu_obj       = seu_sub,
    assay         = "RNA",
    ident1        = g1,
    ident2        = g2,
    features      = features_rna,
    min_pct       = 0,
    logfc         = 0,
    test_use      = "wilcox",
    min_cells     = min_cells_per_group,
    return_thresh = 1
  )

  ## ---- GSEA
  gsea_reactome <- run_fgsea_from_rank(deg_rna_rank, pathways_REACTOME, nperm = 10000, nproc = 1)
  gsea_kegg     <- run_fgsea_from_rank(deg_rna_rank, pathways_KEGG,     nperm = 10000, nproc = 1)
  gsea_gobp     <- run_fgsea_from_rank(deg_rna_rank, pathways_GOBP,     nperm = 10000, nproc = 1)
  gsea_hallmark <- run_fgsea_from_rank(deg_rna_rank, pathways_HALLMARK, nperm = 10000, nproc = 1)

  ## ---- Regulon AUC DE
  Seurat::DefaultAssay(seu_sub) <- "AUC"
  deg_auc <- run_findmarkers(
    seu_obj       = seu_sub,
    assay         = "AUC",
    ident1        = g1,
    ident2        = g2,
    features      = rownames(seu_sub[["AUC"]]),
    min_pct       = 0.1,
    logfc         = 0,
    test_use      = "wilcox",
    min_cells     = min_cells_per_group,
    return_thresh = 1
  )

  ## ---- Regulon Bin DE
  Seurat::DefaultAssay(seu_sub) <- "Bin"
  deg_bin <- run_findmarkers(
    seu_obj       = seu_sub,
    assay         = "Bin",
    ident1        = g1,
    ident2        = g2,
    features      = rownames(seu_sub[["Bin"]]),
    min_pct       = 0.1,
    logfc         = 0,
    test_use      = "wilcox",
    min_cells     = min_cells_per_group,
    return_thresh = 1
  )

  ## ---- Scores (UCell etc.) cell-level
  score_cell <- score_compare_celllevel(
    seu_obj   = seu_sub,
    score_cols= score_cols,
    ident1    = g1,
    ident2    = g2,
    group_col = "group",
    min_cells = min_cells_per_group
  )

  # save per-comparison excel
  sheets <- list(
    Info = tibble::tibble(
      comparison = cmp_nm,
      ident1 = g1,
      ident2 = g2,
      scope = scope_name,
      n_cells_total = ncol(seu_sub),
      n_cells_ident1 = n1,
      n_cells_ident2 = n2
    ),
    QC_cells_by_group  = qc_cells_group,
    QC_cells_by_stage  = qc_cells_stage,
    QC_cells_by_sample = qc_cells_sample,
    DEG_RNA            = deg_rna,
    RNA_rank_for_GSEA  = deg_rna_rank,
    GSEA_REACTOME      = gsea_reactome,
    GSEA_KEGG          = gsea_kegg,
    GSEA_GO_BP         = gsea_gobp,
    GSEA_Hallmark      = gsea_hallmark,
    DEG_Regulon_AUC    = deg_auc,
    DEG_Regulon_Bin    = deg_bin,
    Scores_cellLevel   = score_cell
  )

  out_xlsx <- file.path(out_dir, paste0("Step4_", cmp_nm, "_", scope_name, ".xlsx"))
  write_excel_openxlsx2(sheets, out_xlsx)

  # store for reversal
  res_all[[cmp_nm]] <- list(
    ident1 = g1, ident2 = g2,
    deg_rna = deg_rna,
    deg_auc = deg_auc,
    deg_bin = deg_bin,
    scores  = score_cell,
    gsea = list(
      REACTOME = gsea_reactome,
      KEGG     = gsea_kegg,
      GO_BP    = gsea_gobp,
      HALLMARK = gsea_hallmark
    )
  )

  rm(seu_sub, deg_rna, deg_rna_rank, deg_auc, deg_bin, score_cell,
     gsea_reactome, gsea_kegg, gsea_gobp, gsea_hallmark)
  gc()
}

message("All comparisons finished. Results in: ", out_dir)


## -------------------------- Section 7: Reversal highlights ---------------
# Goal:
#   Base = LEN_vs_CTR  (LEN up => +logFC)
#   Reversal if:
#     - feature significant in base AND significant in (LENnAMG_vs_LEN) with opposite sign
#     - OR significant in base AND significant in (LENnKO_vs_LEN) with opposite sign
#
# NOTE: This is still cell-level pseudo-replication (1 sample per group),
#       but it's useful for candidate screening.

padj_th  <- 0.05

rna_fc_th  <- 0.25
auc_fc_th  <- 0.05
bin_fc_th  <- 0.05
score_th   <- 0.01   # delta_median threshold

gsea_nes_th <- 1.0

make_reversal_table_fc <- function(base_df, mod_df,
                                   base_name, mod_name,
                                   fc_th, padj_th = 0.05) {
  if (is.null(base_df) || is.null(mod_df)) return(NULL)

  fc1 <- detect_fc_col(base_df)
  fc2 <- detect_fc_col(mod_df)

  b <- base_df %>%
    dplyr::filter(!is.na(p_val_adj)) %>%
    dplyr::filter(p_val_adj <= padj_th) %>%
    dplyr::filter(abs(.data[[fc1]]) >= fc_th) %>%
    dplyr::select(feature, p_val_adj, dplyr::all_of(fc1)) %>%
    dplyr::rename(base_padj = p_val_adj, base_fc = dplyr::all_of(fc1))

  m <- mod_df %>%
    dplyr::filter(!is.na(p_val_adj)) %>%
    dplyr::filter(p_val_adj <= padj_th) %>%
    dplyr::filter(abs(.data[[fc2]]) >= fc_th) %>%
    dplyr::select(feature, p_val_adj, dplyr::all_of(fc2)) %>%
    dplyr::rename(mod_padj = p_val_adj, mod_fc = dplyr::all_of(fc2))

  j <- dplyr::inner_join(b, m, by = "feature") %>%
    dplyr::mutate(
      base_cmp = base_name,
      mod_cmp  = mod_name,
      reversed = sign(base_fc) == -sign(mod_fc)
    ) %>%
    dplyr::filter(reversed) %>%
    dplyr::arrange(base_padj, mod_padj, dplyr::desc(abs(base_fc)))

  j
}

make_reversal_table_scores <- function(base_df, mod_df,
                                       base_name, mod_name,
                                       th = 0.01, padj_th = 0.05) {
  if (is.null(base_df) || is.null(mod_df)) return(NULL)

  b <- base_df %>%
    dplyr::filter(!is.na(p_adj)) %>%
    dplyr::filter(p_adj <= padj_th) %>%
    dplyr::filter(abs(delta_median) >= th) %>%
    dplyr::select(score, p_adj, delta_median, median1, median2) %>%
    dplyr::rename(base_padj = p_adj, base_delta = delta_median,
                  base_med1 = median1, base_med2 = median2)

  m <- mod_df %>%
    dplyr::filter(!is.na(p_adj)) %>%
    dplyr::filter(p_adj <= padj_th) %>%
    dplyr::filter(abs(delta_median) >= th) %>%
    dplyr::select(score, p_adj, delta_median, median1, median2) %>%
    dplyr::rename(mod_padj = p_adj, mod_delta = delta_median,
                  mod_med1 = median1, mod_med2 = median2)

  j <- dplyr::inner_join(b, m, by = "score") %>%
    dplyr::mutate(
      base_cmp = base_name,
      mod_cmp  = mod_name,
      reversed = sign(base_delta) == -sign(mod_delta)
    ) %>%
    dplyr::filter(reversed) %>%
    dplyr::arrange(base_padj, mod_padj, dplyr::desc(abs(base_delta)))

  j
}

make_reversal_table_gsea <- function(base_df, mod_df,
                                     base_name, mod_name,
                                     nes_th = 1.0, padj_th = 0.05) {
  if (is.null(base_df) || is.null(mod_df)) return(NULL)
  if (!all(c("pathway","NES","padj") %in% colnames(base_df))) return(NULL)
  if (!all(c("pathway","NES","padj") %in% colnames(mod_df))) return(NULL)

  b <- base_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::filter(padj <= padj_th) %>%
    dplyr::filter(abs(NES) >= nes_th) %>%
    dplyr::select(pathway, NES, padj) %>%
    dplyr::rename(base_NES = NES, base_padj = padj)

  m <- mod_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::filter(padj <= padj_th) %>%
    dplyr::filter(abs(NES) >= nes_th) %>%
    dplyr::select(pathway, NES, padj) %>%
    dplyr::rename(mod_NES = NES, mod_padj = padj)

  j <- dplyr::inner_join(b, m, by = "pathway") %>%
    dplyr::mutate(
      base_cmp = base_name,
      mod_cmp  = mod_name,
      reversed = sign(base_NES) == -sign(mod_NES)
    ) %>%
    dplyr::filter(reversed) %>%
    dplyr::arrange(base_padj, mod_padj, dplyr::desc(abs(base_NES)))

  j
}

# Require base comparison
if (!("LEN_vs_CTR" %in% names(res_all))) {
  warning("LEN_vs_CTR not found in res_all. Reversal summary will be skipped.")
} else {

  base <- res_all[["LEN_vs_CTR"]]

  rev_AMG_rna <- make_reversal_table_fc(base$deg_rna, res_all[["LENnAMG_vs_LEN"]]$deg_rna,
                                       "LEN_vs_CTR", "LENnAMG_vs_LEN", fc_th = rna_fc_th, padj_th = padj_th)
  rev_KO_rna  <- make_reversal_table_fc(base$deg_rna, res_all[["LENnKO_vs_LEN"]]$deg_rna,
                                       "LEN_vs_CTR", "LENnKO_vs_LEN", fc_th = rna_fc_th, padj_th = padj_th)

  rev_AMG_auc <- make_reversal_table_fc(base$deg_auc, res_all[["LENnAMG_vs_LEN"]]$deg_auc,
                                       "LEN_vs_CTR", "LENnAMG_vs_LEN", fc_th = auc_fc_th, padj_th = padj_th)
  rev_KO_auc  <- make_reversal_table_fc(base$deg_auc, res_all[["LENnKO_vs_LEN"]]$deg_auc,
                                       "LEN_vs_CTR", "LENnKO_vs_LEN", fc_th = auc_fc_th, padj_th = padj_th)

  rev_AMG_bin <- make_reversal_table_fc(base$deg_bin, res_all[["LENnAMG_vs_LEN"]]$deg_bin,
                                       "LEN_vs_CTR", "LENnAMG_vs_LEN", fc_th = bin_fc_th, padj_th = padj_th)
  rev_KO_bin  <- make_reversal_table_fc(base$deg_bin, res_all[["LENnKO_vs_LEN"]]$deg_bin,
                                       "LEN_vs_CTR", "LENnKO_vs_LEN", fc_th = bin_fc_th, padj_th = padj_th)

  rev_AMG_score <- make_reversal_table_scores(base$scores, res_all[["LENnAMG_vs_LEN"]]$scores,
                                             "LEN_vs_CTR", "LENnAMG_vs_LEN", th = score_th, padj_th = padj_th)
  rev_KO_score  <- make_reversal_table_scores(base$scores, res_all[["LENnKO_vs_LEN"]]$scores,
                                             "LEN_vs_CTR", "LENnKO_vs_LEN", th = score_th, padj_th = padj_th)

  # GSEA reversal per collection
  gsea_types <- c("REACTOME","KEGG","GO_BP","HALLMARK")
  rev_gsea_AMG <- list()
  rev_gsea_KO  <- list()

  for (tp in gsea_types) {
    rev_gsea_AMG[[tp]] <- make_reversal_table_gsea(
      base$gsea[[tp]],
      res_all[["LENnAMG_vs_LEN"]]$gsea[[tp]],
      "LEN_vs_CTR", "LENnAMG_vs_LEN",
      nes_th = gsea_nes_th, padj_th = padj_th
    )
    rev_gsea_KO[[tp]] <- make_reversal_table_gsea(
      base$gsea[[tp]],
      res_all[["LENnKO_vs_LEN"]]$gsea[[tp]],
      "LEN_vs_CTR", "LENnKO_vs_LEN",
      nes_th = gsea_nes_th, padj_th = padj_th
    )
  }

  # Summary counts
  summary_counts <- tibble::tibble(
    item = c("RNA_genes", "Regulon_AUC", "Regulon_Bin", "Scores",
             paste0("GSEA_", gsea_types)),
    reversed_by_AMG = c(
      ifelse(is.null(rev_AMG_rna), 0, nrow(rev_AMG_rna)),
      ifelse(is.null(rev_AMG_auc), 0, nrow(rev_AMG_auc)),
      ifelse(is.null(rev_AMG_bin), 0, nrow(rev_AMG_bin)),
      ifelse(is.null(rev_AMG_score), 0, nrow(rev_AMG_score)),
      vapply(gsea_types, function(tp) ifelse(is.null(rev_gsea_AMG[[tp]]), 0, nrow(rev_gsea_AMG[[tp]])), integer(1))
    ),
    reversed_by_KO = c(
      ifelse(is.null(rev_KO_rna), 0, nrow(rev_KO_rna)),
      ifelse(is.null(rev_KO_auc), 0, nrow(rev_KO_auc)),
      ifelse(is.null(rev_KO_bin), 0, nrow(rev_KO_bin)),
      ifelse(is.null(rev_KO_score), 0, nrow(rev_KO_score)),
      vapply(gsea_types, function(tp) ifelse(is.null(rev_gsea_KO[[tp]]), 0, nrow(rev_gsea_KO[[tp]])), integer(1))
    )
  )

  # write reversal excel
  sheets_rev <- list(
    Summary_counts = summary_counts,
    RNA_reversed_by_AMG = rev_AMG_rna,
    RNA_reversed_by_KO  = rev_KO_rna,
    AUC_reversed_by_AMG = rev_AMG_auc,
    AUC_reversed_by_KO  = rev_KO_auc,
    Bin_reversed_by_AMG = rev_AMG_bin,
    Bin_reversed_by_KO  = rev_KO_bin,
    Scores_reversed_by_AMG = rev_AMG_score,
    Scores_reversed_by_KO  = rev_KO_score
  )

  # add gsea sheets
  for (tp in gsea_types) {
    sheets_rev[[paste0("GSEA_", tp, "_rev_AMG")]] <- rev_gsea_AMG[[tp]]
    sheets_rev[[paste0("GSEA_", tp, "_rev_KO")]]  <- rev_gsea_KO[[tp]]
  }

  out_rev <- file.path(out_dir, "Step4_ReversalHighlights_AllNeu.xlsx")
  write_excel_openxlsx2(sheets_rev, out_rev)
}

message("Step4 DONE.")
