#!/usr/bin/env Rscript

# Six-condition, lineage-resolved scRNA-seq analysis across eight immune lineages.
#
# Inference limitation:
# Each condition is represented by one pooled 10x library. Cell-level Wilcoxon
# tests and BH-adjusted P values are exploratory distributional summaries, not
# treatment-level biological-replicate inference. Matched-background
# attenuation is an effect-size criterion, not a formal interaction/rescue test.
#
# FindMarkers() uses an explicit LogNormalize mean function, and each fold change
# is independently recalculated from the RNA data matrix for numerical QA.

rm(list = ls())
gc()

script_version <- "1.0.0"
canonical_md5 <- "e0c009949ae4392a3779ec1217dbfbb3"
default_input <- NA_character_
default_outdir <- file.path("results", "Step7_Len_CD45_lineage_remodeling")

usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript R/Analysis/Step7_Len_CD45_lineage_remodeling.R [options]",
    "",
    "Options:",
    "  --input PATH          Canonical CD45+ Seurat RDS (required)",
    "  --outdir PATH         New output directory",
    "  --expected-md5 VALUE  Expected input MD5, or 'none' to disable",
    "  --overwrite           Explicitly allow replacement of named outputs",
    "  --help                Print this help and exit",
    "",
    "Both --name=value and --name value forms are accepted.",
    sep = "\n"
  ))
}

parse_cli <- function(argv) {
  config <- list(
    input = default_input,
    outdir = default_outdir,
    expected_md5 = canonical_md5,
    overwrite = FALSE
  )
  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (identical(arg, "--help")) {
      usage()
      quit(save = "no", status = 0L)
    } else if (identical(arg, "--overwrite")) {
      config$overwrite <- TRUE
    } else if (arg %in% c("--input", "--outdir", "--expected-md5")) {
      if (i == length(argv)) stop("Missing value after ", arg, call. = FALSE)
      i <- i + 1L
      value <- argv[[i]]
      name <- sub("^--", "", arg)
      name <- gsub("-", "_", name, fixed = TRUE)
      config[[name]] <- value
    } else if (grepl("^--input=", arg)) {
      config$input <- sub("^--input=", "", arg)
    } else if (grepl("^--outdir=", arg)) {
      config$outdir <- sub("^--outdir=", "", arg)
    } else if (grepl("^--expected-md5=", arg)) {
      config$expected_md5 <- sub("^--expected-md5=", "", arg)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
    i <- i + 1L
  }
  if (is.na(config$input) || !nzchar(config$input)) {
    stop("--input is required.", call. = FALSE)
  }
  if (!nzchar(config$outdir)) stop("--outdir cannot be empty.", call. = FALSE)
  if (tolower(config$expected_md5) %in% c("none", "na", "off")) {
    config$expected_md5 <- NA_character_
  }
  config
}

config <- parse_cli(commandArgs(trailingOnly = TRUE))

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble",
  "ggplot2", "ggVennDiagram", "patchwork", "clusterProfiler",
  "org.Mm.eg.db", "AnnotationDbi", "GOSemSim", "ragg"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    ". Install them in the analysis environment before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggVennDiagram)
  library(patchwork)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(GOSemSim)
  library(ragg)
})

group_levels <- c("CTR", "LEN", "CTRnAMG", "LENnAMG", "CTRnKO", "LENnKO")
lineage_levels <- c("Neu", "Mac", "Mono", "DC", "CD8T", "NonCD8_T", "NK", "B")
myeloid_levels <- c("Neu", "Mac", "Mono", "DC")
lymphoid_levels <- c("CD8T", "NonCD8_T", "NK", "B")
cd8_l2 <- c(
  "c1_T_CD8eff_Cxcr3", "c2_T_CD8exh_Pdcd1Lag3", "c5_T_CD8cycle_Mki67"
)
non_cd8_t_l2 <- c(
  "c3_T_Treg_Foxp3_Il2ra", "c4_T_Naive_Tcf7_Il7r", "c6_T_NKT_Zbtb16_Tcf7"
)
contrast_definitions <- list(
  WT = c(treated = "LEN", control = "CTR"),
  AMG = c(treated = "LENnAMG", control = "CTRnAMG"),
  KO = c(treated = "LENnKO", control = "CTRnKO")
)

min_cells_per_group <- 30L
min_detection_fraction <- 0.10
deg_abs_log2fc <- 0.25
deg_bh_fdr <- 0.05
attenuation_abs_did <- 0.25
attenuation_fraction <- 0.50
minimum_mappable_center_genes <- 10L

dirs <- list(
  root = config$outdir,
  tables = file.path(config$outdir, "tables"),
  figures = file.path(config$outdir, "figures"),
  logs = file.path(config$outdir, "logs")
)
for (path in unname(unlist(dirs))) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create output directory: ", path, call. = FALSE)
  }
}

output_paths <- list(
  contract = file.path(dirs$root, "Step7_analysis_contract.csv"),
  input_manifest = file.path(dirs$root, "Step7_input_manifest.csv"),
  cell_counts = file.path(dirs$tables, "Step7_lineage_cell_counts_by_group.csv"),
  universe = file.path(dirs$tables, "Step7_lineage_gene_universe.csv"),
  de_long = file.path(dirs$tables, "Step7_lineage_DE_all_contrasts.csv"),
  de_wide = file.path(dirs$tables, "Step7_lineage_DE_wide_effects.csv"),
  fc_qa = file.path(dirs$tables, "Step7_explicit_LogNormalize_FC_QA.csv"),
  classification = file.path(dirs$tables, "Step7_attenuation_classification.csv"),
  membership = file.path(dirs$tables, "Step7_Venn_set_membership.csv"),
  venn_counts = file.path(dirs$tables, "Step7_Venn_counts.csv"),
  center_genes = file.path(dirs$tables, "Step7_joint_center_genes.csv"),
  go_status = file.path(dirs$tables, "Step7_GO_BP_status.csv"),
  go_mapping = file.path(dirs$tables, "Step7_GO_BP_gene_mapping.csv"),
  go_full = file.path(dirs$tables, "Step7_GO_BP_full.csv"),
  go_simplified = file.path(dirs$tables, "Step7_GO_BP_simplified_significant.csv"),
  go_display = file.path(dirs$tables, "Step7_GO_BP_display_terms.csv"),
  qa = file.path(dirs$root, "Step7_QA_checks.csv"),
  report = file.path(dirs$root, "Step7_summary.md"),
  session = file.path(dirs$logs, "Step7_sessionInfo.txt"),
  fig13_png = file.path(dirs$figures, "Figure13_myeloid_Venn_GO.png"),
  fig13_pdf = file.path(dirs$figures, "Figure13_myeloid_Venn_GO.pdf"),
  fig17_png = file.path(dirs$figures, "Figure17_lymphoid_Venn_GO.png"),
  fig17_pdf = file.path(dirs$figures, "Figure17_lymphoid_Venn_GO.pdf"),
  output_manifest = file.path(dirs$root, "Step7_output_manifest.csv")
)

existing_targets <- unlist(output_paths)[file.exists(unlist(output_paths))]
if (length(existing_targets) && !config$overwrite) {
  stop(
    "Refusing to overwrite existing named output(s). Re-run with --overwrite only after review:\n",
    paste(existing_targets, collapse = "\n"), call. = FALSE
  )
}

publish_temp <- function(temp_path, final_path) {
  if (!file.exists(temp_path)) stop("Temporary output was not created: ", temp_path, call. = FALSE)
  if (file.exists(final_path) && !config$overwrite) {
    stop("Refusing to overwrite: ", final_path, call. = FALSE)
  }
  backup_path <- NA_character_
  if (file.exists(final_path)) {
    backup_path <- paste0(final_path, ".backup-", Sys.getpid())
    if (file.exists(backup_path)) stop("Backup path already exists: ", backup_path, call. = FALSE)
    if (!file.rename(final_path, backup_path)) {
      stop("Could not stage existing output for explicit overwrite: ", final_path, call. = FALSE)
    }
  }
  if (!file.rename(temp_path, final_path)) {
    if (!is.na(backup_path) && file.exists(backup_path)) file.rename(backup_path, final_path)
    stop("Atomic publish failed: ", final_path, call. = FALSE)
  }
  if (!is.na(backup_path) && file.exists(backup_path)) unlink(backup_path)
  normalizePath(final_path, winslash = "/", mustWork = TRUE)
}

write_csv_atomic <- function(x, path) {
  temp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)
  utils::write.csv(x, temp, row.names = FALSE, na = "", quote = TRUE)
  publish_temp(temp, path)
}

write_lines_atomic <- function(x, path) {
  temp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)
  writeLines(x, temp, useBytes = TRUE)
  publish_temp(temp, path)
}

save_plot_atomic <- function(plot, path, width, height, dpi = 320) {
  temp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)
  extension <- tolower(tools::file_ext(path))
  if (identical(extension, "png")) {
    ggplot2::ggsave(
      temp, plot = plot, width = width, height = height, units = "in",
      dpi = dpi, bg = "white", device = ragg::agg_png, limitsize = FALSE
    )
  } else if (identical(extension, "pdf")) {
    ggplot2::ggsave(
      temp, plot = plot, width = width, height = height, units = "in",
      bg = "white", device = grDevices::cairo_pdf, limitsize = FALSE
    )
  } else {
    stop("Unsupported figure extension: ", extension, call. = FALSE)
  }
  publish_temp(temp, path)
}

join_assay_layers_if_needed <- function(object, assay = "RNA") {
  if (!assay %in% names(object@assays)) stop("Missing assay: ", assay, call. = FALSE)
  assay_object <- object[[assay]]
  if (!inherits(assay_object, "Assay5")) return(object)
  layers <- SeuratObject::Layers(assay_object)
  layer_families <- sub("[.].*$", "", layers)
  if (anyDuplicated(layer_families)) {
    object <- SeuratObject::JoinLayers(object, assay = assay)
  }
  object
}

get_assay_matrix <- function(object, assay = "RNA", layer = "counts") {
  if (inherits(object[[assay]], "Assay5")) {
    layers <- SeuratObject::Layers(object[[assay]])
    if (!layer %in% layers) stop("Missing RNA layer after JoinLayers: ", layer, call. = FALSE)
    return(SeuratObject::LayerData(object, assay = assay, layer = layer))
  }
  Seurat::GetAssayData(object, assay = assay, slot = layer)
}

map_symbols <- function(symbols) {
  symbols <- unique(stats::na.omit(as.character(symbols)))
  symbols <- symbols[nzchar(symbols)]
  if (!length(symbols)) {
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  }
  suppressMessages(
    AnnotationDbi::select(
      org.Mm.eg.db::org.Mm.eg.db,
      keys = symbols,
      keytype = "SYMBOL",
      columns = "ENTREZID"
    )
  ) |>
    dplyr::filter(!is.na(ENTREZID)) |>
    dplyr::distinct(SYMBOL, ENTREZID)
}

parse_gene_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(parts) {
    if (length(parts) != 2L) return(NA_real_)
    as.numeric(parts[[1]]) / as.numeric(parts[[2]])
  }, numeric(1))
}

theme_clean <- function(base_size = 8) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "#777777", linewidth = 0.4),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      axis.text = ggplot2::element_text(color = "#333333"),
      plot.margin = ggplot2::margin(5, 7, 5, 5)
    )
}

if (!file.exists(config$input) || isTRUE(file.info(config$input)$isdir)) {
  stop("Input Seurat RDS is missing: ", config$input, call. = FALSE)
}
input_md5 <- unname(tools::md5sum(config$input))
if (!is.na(config$expected_md5) && !identical(tolower(input_md5), tolower(config$expected_md5))) {
  stop(
    "Input MD5 mismatch. Observed ", input_md5, "; expected ", config$expected_md5,
    ". Use --expected-md5 none only for an intentionally changed input.", call. = FALSE
  )
}

message("Reading canonical CD45+ Seurat object: ", config$input)
seu <- readRDS(config$input)
if (!inherits(seu, "Seurat")) stop("Input is not a Seurat object.", call. = FALSE)
if (!"RNA" %in% SeuratObject::Assays(seu)) stop("Input lacks the RNA assay.", call. = FALSE)
seu <- join_assay_layers_if_needed(seu, "RNA")
SeuratObject::DefaultAssay(seu) <- "RNA"
counts <- get_assay_matrix(seu, "RNA", "counts")
data_matrix <- get_assay_matrix(seu, "RNA", "data")
if (!identical(dim(counts), dim(data_matrix)) ||
    !identical(rownames(counts), rownames(data_matrix)) ||
    !identical(colnames(counts), colnames(data_matrix))) {
  stop("RNA counts and data layers are not dimensionally identical.", call. = FALSE)
}

required_meta <- c("group", "rev_CellType_l1", "rev_CellType_l2")
missing_meta <- setdiff(required_meta, colnames(seu[[]]))
if (length(missing_meta)) {
  stop("Missing required metadata: ", paste(missing_meta, collapse = ", "), call. = FALSE)
}

meta <- seu[[]] |>
  tibble::rownames_to_column("cell") |>
  dplyr::transmute(
    cell,
    group = as.character(group),
    rev_CellType_l1 = as.character(rev_CellType_l1),
    rev_CellType_l2 = as.character(rev_CellType_l2),
    lineage = dplyr::case_when(
      rev_CellType_l2 %in% cd8_l2 ~ "CD8T",
      rev_CellType_l2 %in% non_cd8_t_l2 ~ "NonCD8_T",
      rev_CellType_l1 %in% c("Neu", "Mac", "Mono", "DC", "NK", "B") ~ rev_CellType_l1,
      TRUE ~ NA_character_
    )
  )
if (!identical(meta$cell, colnames(counts))) {
  meta <- meta[match(colnames(counts), meta$cell), , drop = FALSE]
}
if (anyNA(meta$cell) || !identical(meta$cell, colnames(counts))) {
  stop("Metadata and RNA matrix cell identifiers do not match.", call. = FALSE)
}
if (!setequal(unique(meta$group), group_levels)) {
  stop(
    "Unexpected group labels. Observed: ", paste(sort(unique(meta$group)), collapse = ", "),
    call. = FALSE
  )
}
unknown_t <- unique(meta$rev_CellType_l2[
  meta$rev_CellType_l1 == "T" &
    !meta$rev_CellType_l2 %in% c(cd8_l2, non_cd8_t_l2)
])
unknown_t <- stats::na.omit(unknown_t)
if (length(unknown_t)) {
  stop("Unmapped frozen T-cell labels: ", paste(unknown_t, collapse = ", "), call. = FALSE)
}

meta$group <- factor(meta$group, levels = group_levels)
cell_counts <- meta |>
  dplyr::filter(!is.na(lineage)) |>
  dplyr::count(lineage, group, name = "n_cells") |>
  tidyr::complete(
    lineage = lineage_levels,
    group = group_levels,
    fill = list(n_cells = 0L)
  ) |>
  dplyr::mutate(
    lineage = factor(lineage, levels = lineage_levels),
    group = factor(group, levels = group_levels)
  ) |>
  dplyr::arrange(lineage, group) |>
  dplyr::mutate(
    eligible_ge_30 = n_cells >= min_cells_per_group,
    lineage = as.character(lineage),
    group = as.character(group)
  )
if (!all(cell_counts$eligible_ge_30)) {
  failed <- cell_counts |>
    dplyr::filter(!eligible_ge_30) |>
    dplyr::transmute(key = paste(lineage, group, n_cells, sep = ":")) |>
    dplyr::pull(key)
  stop("Fewer than 30 cells in required lineage/group: ", paste(failed, collapse = ", "), call. = FALSE)
}

meta$analysis_identity <- ifelse(
  is.na(meta$lineage),
  paste0("Excluded__", as.character(meta$group)),
  paste0(meta$lineage, "__", as.character(meta$group))
)
identity_by_cell <- stats::setNames(meta$analysis_identity, meta$cell)
seu$Step7_analysis_identity <- unname(identity_by_cell[colnames(seu)])
SeuratObject::Idents(seu) <- "Step7_analysis_identity"
expected_identities <- as.vector(outer(lineage_levels, group_levels, paste, sep = "__"))
missing_identities <- setdiff(expected_identities, levels(SeuratObject::Idents(seu)))
if (length(missing_identities)) {
  stop("Missing lineage/group identities: ", paste(missing_identities, collapse = ", "), call. = FALSE)
}

lognormalize_fc_mean <- function(x) {
  if (!ncol(x)) stop("Cannot calculate a group mean from zero cells.", call. = FALSE)
  log2((Matrix::rowSums(expm1(x)) + 1) / ncol(x))
}

fc_qa_rows <- list()

find_one_contrast <- function(lineage_name, eligible_genes, contrast_name, definition) {
  treated_ident <- paste0(lineage_name, "__", definition[["treated"]])
  control_ident <- paste0(lineage_name, "__", definition[["control"]])
  n_treated <- sum(meta$lineage == lineage_name & meta$group == definition[["treated"]], na.rm = TRUE)
  n_control <- sum(meta$lineage == lineage_name & meta$group == definition[["control"]], na.rm = TRUE)
  result <- Seurat::FindMarkers(
    object = seu,
    ident.1 = treated_ident,
    ident.2 = control_ident,
    assay = "RNA",
    features = eligible_genes,
    test.use = "wilcox",
    slot = "data",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE,
    mean.fxn = lognormalize_fc_mean,
    fc.name = "avg_log2FC",
    base = 2,
    densify = FALSE,
    verbose = FALSE
  ) |>
    tibble::rownames_to_column("gene")
  fc_candidates <- c("avg_log2FC", "avg_logFC", "avg_diff")
  fc_column <- fc_candidates[fc_candidates %in% colnames(result)]
  if (!length(fc_column)) stop("Cannot identify a log-fold-change column.", call. = FALSE)
  if (!identical(fc_column[[1]], "avg_log2FC")) {
    result <- dplyr::rename(result, avg_log2FC = dplyr::all_of(fc_column[[1]]))
  }
  result <- result |>
    dplyr::mutate(
      lineage = lineage_name,
      contrast = contrast_name,
      treated = unname(definition[["treated"]]),
      control = unname(definition[["control"]]),
      n_treated = n_treated,
      n_control = n_control,
      p_adj_BH = stats::p.adjust(p_val, method = "BH"),
      direction = dplyr::case_when(
        avg_log2FC > deg_abs_log2fc & p_adj_BH < deg_bh_fdr ~ "Up",
        avg_log2FC < -deg_abs_log2fc & p_adj_BH < deg_bh_fdr ~ "Down",
        TRUE ~ "Not DEG"
      )
    ) |>
    dplyr::select(
      lineage, contrast, treated, control, gene, avg_log2FC,
      p_val, p_adj_BH, pct.1, pct.2, n_treated, n_control, direction,
      dplyr::everything()
    )
  if (!setequal(result$gene, eligible_genes)) {
    stop(
      "FindMarkers did not return the complete common universe for ",
      lineage_name, " / ", contrast_name, call. = FALSE
    )
  }
  treated_cells <- meta$cell[which(
    meta$lineage == lineage_name & meta$group == definition[["treated"]]
  )]
  control_cells <- meta$cell[which(
    meta$lineage == lineage_name & meta$group == definition[["control"]]
  )]
  manual_fc <- lognormalize_fc_mean(
    data_matrix[eligible_genes, treated_cells, drop = FALSE]
  ) - lognormalize_fc_mean(
    data_matrix[eligible_genes, control_cells, drop = FALSE]
  )
  names(manual_fc) <- eligible_genes
  fc_error <- result$avg_log2FC - unname(manual_fc[result$gene])
  qa_key <- paste(lineage_name, contrast_name, sep = "__")
  fc_qa_rows[[qa_key]] <<- data.frame(
    lineage = lineage_name,
    contrast = contrast_name,
    n_genes = length(eligible_genes),
    max_abs_formula_error = max(abs(fc_error), na.rm = TRUE),
    median_abs_formula_error = stats::median(abs(fc_error), na.rm = TRUE),
    Ltf_FindMarkers_log2FC = if ("Ltf" %in% result$gene) {
      result$avg_log2FC[result$gene == "Ltf"][[1]]
    } else NA_real_,
    Ltf_manual_log2FC = if ("Ltf" %in% names(manual_fc)) {
      manual_fc[["Ltf"]]
    } else NA_real_,
    stringsAsFactors = FALSE
  )
  result
}

universe_rows <- list()
de_rows <- list()
for (lineage_name in lineage_levels) {
  message("Differential expression: ", lineage_name)
  group_indices <- lapply(group_levels, function(group_name) {
    which(meta$lineage == lineage_name & meta$group == group_name)
  })
  names(group_indices) <- group_levels
  detection <- vapply(group_levels, function(group_name) {
    Matrix::rowMeans(counts[, group_indices[[group_name]], drop = FALSE] > 0)
  }, numeric(nrow(counts)))
  rownames(detection) <- rownames(counts)
  eligible <- rownames(counts)[apply(detection >= min_detection_fraction, 1L, any)]
  if (!length(eligible)) stop("No eligible genes for lineage: ", lineage_name, call. = FALSE)
  detection_table <- as.data.frame(detection[eligible, , drop = FALSE], check.names = FALSE) |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(
      lineage = lineage_name,
      max_detection_fraction = apply(detection[eligible, , drop = FALSE], 1L, max),
      .before = 1L
    )
  universe_rows[[lineage_name]] <- detection_table
  for (contrast_name in names(contrast_definitions)) {
    key <- paste(lineage_name, contrast_name, sep = "__")
    de_rows[[key]] <- find_one_contrast(
      lineage_name,
      eligible,
      contrast_name,
      contrast_definitions[[contrast_name]]
    )
  }
}

gene_universe <- dplyr::bind_rows(universe_rows) |>
  dplyr::arrange(match(lineage, lineage_levels), gene)
de_all <- dplyr::bind_rows(de_rows) |>
  dplyr::arrange(match(lineage, lineage_levels), match(contrast, names(contrast_definitions)), gene)
fc_qa <- dplyr::bind_rows(fc_qa_rows) |>
  dplyr::arrange(match(lineage, lineage_levels), match(contrast, names(contrast_definitions)))
de_wide <- de_all |>
  dplyr::select(
    lineage, gene, contrast, avg_log2FC, p_val, p_adj_BH, pct.1, pct.2, direction
  ) |>
  tidyr::pivot_wider(
    names_from = contrast,
    values_from = c(avg_log2FC, p_val, p_adj_BH, pct.1, pct.2, direction),
    names_glue = "{contrast}_{.value}"
  ) |>
  dplyr::arrange(match(lineage, lineage_levels), gene)

effect_table <- de_wide |>
  dplyr::mutate(
    AMG_DiD = AMG_avg_log2FC - WT_avg_log2FC,
    KO_DiD = KO_avg_log2FC - WT_avg_log2FC,
    WT_direction_sign = sign(WT_avg_log2FC),
    AMG_opposing_DiD = WT_direction_sign * AMG_DiD < 0 &
      abs(AMG_DiD) > attenuation_abs_did,
    KO_opposing_DiD = WT_direction_sign * KO_DiD < 0 &
      abs(KO_DiD) > attenuation_abs_did,
    AMG_direction_flip = (WT_avg_log2FC > 0 & AMG_avg_log2FC <= 0) |
      (WT_avg_log2FC < 0 & AMG_avg_log2FC >= 0),
    KO_direction_flip = (WT_avg_log2FC > 0 & KO_avg_log2FC <= 0) |
      (WT_avg_log2FC < 0 & KO_avg_log2FC >= 0),
    AMG_reduced_at_least_50pct = abs(AMG_avg_log2FC) <=
      attenuation_fraction * abs(WT_avg_log2FC),
    KO_reduced_at_least_50pct = abs(KO_avg_log2FC) <=
      attenuation_fraction * abs(WT_avg_log2FC),
    AMG_attenuated = WT_direction_sign != 0 & AMG_opposing_DiD &
      (AMG_direction_flip | AMG_reduced_at_least_50pct),
    KO_attenuated = WT_direction_sign != 0 & KO_opposing_DiD &
      (KO_direction_flip | KO_reduced_at_least_50pct),
    LEN_up_DEG = WT_avg_log2FC > deg_abs_log2fc & WT_p_adj_BH < deg_bh_fdr,
    LEN_down_DEG = WT_avg_log2FC < -deg_abs_log2fc & WT_p_adj_BH < deg_bh_fdr,
    AMG_attenuates_LEN_up = WT_avg_log2FC > 0 & AMG_attenuated,
    KO_attenuates_LEN_up = WT_avg_log2FC > 0 & KO_attenuated,
    AMG_attenuates_LEN_down = WT_avg_log2FC < 0 & AMG_attenuated,
    KO_attenuates_LEN_down = WT_avg_log2FC < 0 & KO_attenuated,
    joint_center_up = LEN_up_DEG & AMG_attenuates_LEN_up & KO_attenuates_LEN_up,
    joint_center_down = LEN_down_DEG & AMG_attenuates_LEN_down & KO_attenuates_LEN_down
  )

membership <- dplyr::bind_rows(
  effect_table |>
    dplyr::transmute(
      lineage, direction = "Up", gene,
      LEN_DEG = LEN_up_DEG,
      AMG_attenuation = AMG_attenuates_LEN_up,
      KO_attenuation = KO_attenuates_LEN_up,
      joint_center = joint_center_up,
      WT_avg_log2FC, AMG_avg_log2FC, KO_avg_log2FC,
      AMG_DiD, KO_DiD, WT_p_adj_BH
    ),
  effect_table |>
    dplyr::transmute(
      lineage, direction = "Down", gene,
      LEN_DEG = LEN_down_DEG,
      AMG_attenuation = AMG_attenuates_LEN_down,
      KO_attenuation = KO_attenuates_LEN_down,
      joint_center = joint_center_down,
      WT_avg_log2FC, AMG_avg_log2FC, KO_avg_log2FC,
      AMG_DiD, KO_DiD, WT_p_adj_BH
    )
) |>
  dplyr::arrange(match(lineage, lineage_levels), match(direction, c("Up", "Down")), gene)

# Distinct summary names prevent within-summarise name reuse.
venn_counts <- membership |>
  dplyr::group_by(lineage, direction) |>
  dplyr::summarise(
    tested_universe_n = dplyr::n(),
    LEN_DEG_n = sum(LEN_DEG),
    AMG_attenuation_n = sum(AMG_attenuation),
    KO_attenuation_n = sum(KO_attenuation),
    LEN_and_AMG_n = sum(LEN_DEG & AMG_attenuation),
    LEN_and_KO_n = sum(LEN_DEG & KO_attenuation),
    AMG_and_KO_n = sum(AMG_attenuation & KO_attenuation),
    joint_center_n = sum(LEN_DEG & AMG_attenuation & KO_attenuation),
    .groups = "drop"
  ) |>
  dplyr::arrange(match(lineage, lineage_levels), match(direction, c("Up", "Down")))

center_genes <- membership |>
  dplyr::filter(joint_center) |>
  dplyr::arrange(
    match(lineage, lineage_levels), match(direction, c("Up", "Down")),
    WT_p_adj_BH, dplyr::desc(abs(WT_avg_log2FC)), gene
  )

# Prespecified display terms receive priority only when they remain significant.
display_term_manifest <- tibble::tribble(
  ~lineage, ~direction, ~display_order, ~Description,
  "Neu", "Up", 1L, "oxidative phosphorylation",
  "Neu", "Up", 2L, "aerobic respiration",
  "Neu", "Up", 3L, "ATP metabolic process",
  "Neu", "Up", 4L, "proton motive force-driven ATP synthesis",
  "Neu", "Up", 5L, "small molecule metabolic process",
  "Neu", "Up", 6L, "antimicrobial humoral response",
  "Neu", "Down", 1L, "defense response to other organism",
  "Neu", "Down", 2L, "regulation of innate immune response",
  "Neu", "Down", 3L, "response to virus",
  "Neu", "Down", 4L, "interferon-mediated signaling pathway",
  "Neu", "Down", 5L, "antigen processing and presentation",
  "Neu", "Down", 6L, "antigen processing and presentation of endogenous peptide antigen via MHC class I",
  "Mac", "Up", 1L, "positive regulation of cell death",
  "Mac", "Up", 2L, "acute-phase response",
  "Mac", "Up", 3L, "fatty acid transport",
  "Mac", "Down", 1L, "innate immune response",
  "Mac", "Down", 2L, "response to cytokine",
  "Mac", "Down", 3L, "antigen processing and presentation",
  "Mac", "Down", 4L, "immune effector process",
  "Mac", "Down", 5L, "antimicrobial humoral immune response mediated by antimicrobial peptide",
  "Mac", "Down", 6L, "regulation of immune response",
  "Mono", "Up", 1L, "mitochondrion organization",
  "Mono", "Up", 2L, "purine-containing compound metabolic process",
  "Mono", "Up", 3L, "ATP metabolic process",
  "Mono", "Up", 4L, "aerobic respiration",
  "Mono", "Up", 5L, "purine nucleoside triphosphate biosynthetic process",
  "Mono", "Up", 6L, "ribonucleoside triphosphate biosynthetic process",
  "DC", "Up", 1L, "cellular respiration",
  "DC", "Up", 2L, "oxidative phosphorylation",
  "DC", "Up", 3L, "ATP metabolic process",
  "DC", "Up", 4L, "proton motive force-driven ATP synthesis",
  "DC", "Up", 5L, "ribonucleoside triphosphate biosynthetic process",
  "DC", "Up", 6L, "purine nucleoside triphosphate biosynthetic process",
  "DC", "Down", 1L, "regulation of vesicle-mediated transport",
  "CD8T", "Up", 1L, "aerobic respiration",
  "CD8T", "Up", 2L, "oxidative phosphorylation",
  "CD8T", "Up", 3L, "ATP metabolic process",
  "CD8T", "Up", 4L, "proton motive force-driven ATP synthesis",
  "CD8T", "Up", 5L, "nucleoside triphosphate biosynthetic process",
  "CD8T", "Up", 6L, "ribonucleoside triphosphate biosynthetic process",
  "CD8T", "Down", 1L, "inflammatory response",
  "CD8T", "Down", 2L, "cytokine-mediated signaling pathway",
  "CD8T", "Down", 3L, "cytokine production",
  "CD8T", "Down", 4L, "immune effector process",
  "CD8T", "Down", 5L, "leukocyte activation involved in immune response",
  "CD8T", "Down", 6L, "regulation of T-helper 1 type immune response",
  "NonCD8_T", "Up", 1L, "cytoplasmic translation",
  "NonCD8_T", "Up", 2L, "ATP metabolic process",
  "NonCD8_T", "Up", 3L, "cellular respiration",
  "NonCD8_T", "Up", 4L, "proton motive force-driven ATP synthesis",
  "NonCD8_T", "Up", 5L, "purine nucleoside triphosphate biosynthetic process",
  "NonCD8_T", "Up", 6L, "purine ribonucleotide metabolic process",
  "NonCD8_T", "Down", 1L, "cytokine production",
  "NonCD8_T", "Down", 2L, "response to cytokine",
  "NonCD8_T", "Down", 3L, "regulation of cytokine production",
  "NonCD8_T", "Down", 4L, "regulation of T cell proliferation",
  "NonCD8_T", "Down", 5L, "regulation of leukocyte proliferation",
  "NonCD8_T", "Down", 6L, "cell adhesion",
  "NK", "Up", 1L, "aerobic respiration",
  "NK", "Up", 2L, "oxidative phosphorylation",
  "NK", "Up", 3L, "ATP metabolic process",
  "NK", "Up", 4L, "purine nucleoside triphosphate biosynthetic process",
  "NK", "Up", 5L, "ribonucleoside triphosphate biosynthetic process",
  "NK", "Up", 6L, "ribose phosphate biosynthetic process",
  "B", "Up", 1L, "oxidative phosphorylation",
  "B", "Up", 2L, "aerobic respiration",
  "B", "Up", 3L, "proton motive force-driven ATP synthesis",
  "B", "Up", 4L, "mitochondrion organization",
  "B", "Up", 5L, "positive regulation of intrinsic apoptotic signaling pathway",
  "B", "Up", 6L, "intrinsic apoptotic signaling pathway by p53 class mediator",
  "B", "Down", 1L, "glycoprotein metabolic process",
  "B", "Down", 2L, "intracellular protein transport",
  "B", "Down", 3L, "cellular response to growth factor stimulus",
  "B", "Down", 4L, "vesicle organization"
)

message("Preparing mouse GO semantic data.")
sem_data <- GOSemSim::godata("org.Mm.eg.db", ont = "BP", computeIC = FALSE)
go_status_rows <- list()
go_mapping_rows <- list()
go_full_rows <- list()
go_simplified_rows <- list()

for (lineage_name in lineage_levels) {
  universe <- gene_universe |>
    dplyr::filter(lineage == lineage_name) |>
    dplyr::pull(gene) |>
    unique()
  universe_map <- map_symbols(universe)
  for (direction_name in c("Up", "Down")) {
    key <- paste(lineage_name, direction_name, sep = "__")
    genes <- center_genes |>
      dplyr::filter(lineage == lineage_name, direction == direction_name) |>
      dplyr::pull(gene) |>
      unique()
    gene_map <- map_symbols(genes)
    n_mappable_genes <- dplyr::n_distinct(gene_map$SYMBOL)
    gate <- if (n_mappable_genes >= minimum_mappable_center_genes) {
      "GO_ELIGIBLE"
    } else if (n_mappable_genes >= 5L) {
      "LIST_ONLY_5_TO_9_MAPPABLE"
    } else {
      "TOO_FEW_BELOW_5_MAPPABLE"
    }
    enrichment_status <- gate
    full <- data.frame()
    simplified <- data.frame()
    if (identical(gate, "GO_ELIGIBLE")) {
      ego <- suppressMessages(clusterProfiler::enrichGO(
        gene = unique(gene_map$ENTREZID),
        universe = unique(universe_map$ENTREZID),
        OrgDb = org.Mm.eg.db::org.Mm.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
      ))
      full <- as.data.frame(ego)
      if (!nrow(full)) {
        enrichment_status <- "NO_ENRICHMENT_RESULT"
      } else {
        full <- full |>
          dplyr::mutate(GeneRatio_numeric = parse_gene_ratio(GeneRatio)) |>
          dplyr::arrange(p.adjust, dplyr::desc(Count), ID)
        significant_ids <- full |>
          dplyr::filter(is.finite(p.adjust), p.adjust < deg_bh_fdr) |>
          dplyr::pull(ID)
        if (!length(significant_ids)) {
          enrichment_status <- "NO_SIGNIFICANT_GO_BP"
        } else {
          ego_significant <- ego
          ego_significant@result <- ego_significant@result[
            ego_significant@result$ID %in% significant_ids,
            , drop = FALSE
          ]
          simplified_object <- tryCatch(
            clusterProfiler::simplify(
              ego_significant,
              cutoff = 0.70,
              by = "p.adjust",
              select_fun = min,
              measure = "Wang",
              semData = sem_data
            ),
            error = function(e) {
              stop(
                "GO semantic simplify failed for ", key, ": ", conditionMessage(e),
                call. = FALSE
              )
            }
          )
          simplified <- as.data.frame(simplified_object) |>
            dplyr::mutate(GeneRatio_numeric = parse_gene_ratio(GeneRatio)) |>
            dplyr::arrange(p.adjust, dplyr::desc(Count), ID)
          enrichment_status <- "PASS"
        }
      }
    }
    go_status_rows[[key]] <- data.frame(
      lineage = lineage_name,
      direction = direction_name,
      n_center_genes = length(genes),
      n_mappable_center_genes = n_mappable_genes,
      n_universe = length(universe),
      n_mappable_universe_genes = dplyr::n_distinct(universe_map$SYMBOL),
      gate = gate,
      enrichment_status = enrichment_status,
      n_GO_full = nrow(full),
      n_GO_simplified_significant = nrow(simplified),
      stringsAsFactors = FALSE
    )
    if (nrow(gene_map)) {
      go_mapping_rows[[key]] <- gene_map |>
        dplyr::mutate(lineage = lineage_name, direction = direction_name, .before = 1L)
    }
    if (nrow(full)) {
      go_full_rows[[key]] <- full |>
        dplyr::mutate(lineage = lineage_name, direction = direction_name, .before = 1L)
    }
    if (nrow(simplified)) {
      go_simplified_rows[[key]] <- simplified |>
        dplyr::mutate(lineage = lineage_name, direction = direction_name, .before = 1L)
    }
  }
}

go_status <- dplyr::bind_rows(go_status_rows) |>
  dplyr::arrange(match(lineage, lineage_levels), match(direction, c("Up", "Down")))
go_mapping <- dplyr::bind_rows(go_mapping_rows)
go_full <- dplyr::bind_rows(go_full_rows)
go_simplified <- dplyr::bind_rows(go_simplified_rows)

go_display <- go_simplified |>
  dplyr::filter(is.finite(p.adjust), p.adjust < deg_bh_fdr) |>
  dplyr::left_join(
    display_term_manifest |>
      dplyr::rename(prespecified_display_order = display_order),
    by = c("lineage", "direction", "Description"),
    relationship = "one-to-one"
  ) |>
  dplyr::mutate(
    retained_prespecified_priority = !is.na(prespecified_display_order),
    prespecified_display_order = dplyr::coalesce(prespecified_display_order, 9999L)
  ) |>
  dplyr::group_by(lineage, direction) |>
  dplyr::arrange(
    dplyr::desc(retained_prespecified_priority), prespecified_display_order,
    p.adjust, dplyr::desc(Count), ID,
    .by_group = TRUE
  ) |>
  dplyr::slice_head(n = 6L) |>
  dplyr::mutate(
    display_order = dplyr::row_number(),
    display_selection = dplyr::if_else(
      retained_prespecified_priority,
      "prespecified_term_still_significant",
      "statistical_fill_from_simplified_result"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    GeneRatio_numeric = parse_gene_ratio(GeneRatio),
    neg_log10_adjusted_P = -log10(pmax(p.adjust, .Machine$double.xmin)),
    Description_wrapped = vapply(
      Description,
      function(x) paste(strwrap(x, width = 39), collapse = "\n"),
      character(1)
    )
  ) |>
  dplyr::arrange(
    match(lineage, lineage_levels), match(direction, c("Up", "Down")), display_order
  )

make_venn <- function(lineage_name, direction_name) {
  table <- membership |>
    dplyr::filter(lineage == lineage_name, direction == direction_name)
  sets <- list(
    `LEN DEG` = table$gene[table$LEN_DEG],
    `AMG atten.` = table$gene[table$AMG_attenuation],
    `KO atten.` = table$gene[table$KO_attenuation]
  )
  center_n <- sum(table$joint_center)
  high_color <- if (identical(direction_name, "Up")) "#D98C80" else "#7FA4C4"
  suppressWarnings(
    ggVennDiagram::ggVennDiagram(
      sets,
      label = "count",
      label_alpha = 0,
      edge_size = 0.65,
      set_size = 3.1,
      label_size = 3.5
    )
  ) +
    ggplot2::scale_fill_gradient(low = "#FFFFFF", high = high_color) +
    ggplot2::labs(
      title = paste0("LEN-", tolower(direction_name), " genes"),
      subtitle = paste0("Joint center: n=", center_n),
      fill = "Genes"
    ) +
    ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(color = "#444444", hjust = 0),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

color_limits <- range(go_display$neg_log10_adjusted_P, finite = TRUE)
make_go_plot <- function(lineage_name, direction_name) {
  table <- go_display |>
    dplyr::filter(lineage == lineage_name, direction == direction_name) |>
    dplyr::arrange(GeneRatio_numeric, p.adjust) |>
    dplyr::mutate(y_id = factor(seq_len(dplyr::n()), levels = seq_len(dplyr::n())))
  title <- paste0("Joint center GO:BP | LEN-", tolower(direction_name))
  if (!nrow(table)) {
    status <- go_status |>
      dplyr::filter(lineage == lineage_name, direction == direction_name)
    label <- dplyr::case_when(
      status$gate == "LIST_ONLY_5_TO_9_MAPPABLE" ~ "5-9 mappable genes: list only",
      status$gate == "TOO_FEW_BELOW_5_MAPPABLE" ~ "Fewer than 5 mappable genes",
      status$enrichment_status == "NO_SIGNIFICANT_GO_BP" ~ "No significant GO:BP term",
      TRUE ~ "No significant GO:BP display term"
    )
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = label, size = 3, color = "#555555") +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::labs(title = title, x = NULL, y = NULL) +
        ggplot2::theme_void(base_size = 8) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0))
    )
  }
  labels <- stats::setNames(table$Description_wrapped, as.character(table$y_id))
  ggplot2::ggplot(table, ggplot2::aes(GeneRatio_numeric, y_id, fill = neg_log10_adjusted_P)) +
    ggplot2::geom_col(width = 0.72, color = "#666666", linewidth = 0.20) +
    ggplot2::geom_text(
      ggplot2::aes(label = Count), hjust = -0.14, size = 2.35, color = "#222222"
    ) +
    ggplot2::scale_y_discrete(labels = labels) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::scale_fill_viridis_c(
      option = "C", limits = color_limits, name = "-log10 adjusted P"
    ) +
    ggplot2::labs(title = title, x = "GeneRatio", y = NULL) +
    theme_clean(7.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        color = if (identical(direction_name, "Up")) "#B64A3A" else "#356A9A"
      ),
      axis.text.y = ggplot2::element_text(size = 6.25, lineheight = 0.88),
      legend.position = "right"
    )
}

make_lineage_panel <- function(lineage_name, panel_letter) {
  venn_row <- make_venn(lineage_name, "Up") + make_venn(lineage_name, "Down") +
    patchwork::plot_layout(ncol = 2)
  go_row <- make_go_plot(lineage_name, "Up") + make_go_plot(lineage_name, "Down") +
    patchwork::plot_layout(ncol = 2, guides = "collect")
  (venn_row / go_row) +
    patchwork::plot_layout(heights = c(0.92, 1.38)) +
    patchwork::plot_annotation(
      title = paste0(panel_letter, "   ", gsub("_", "-", lineage_name)),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 11.5, hjust = 0)
      )
    ) &
    ggplot2::theme(legend.position = "right")
}

make_figure <- function(lineages, title, subtitle) {
  panels <- lapply(seq_along(lineages), function(index) {
    patchwork::wrap_elements(
      full = make_lineage_panel(lineages[[index]], LETTERS[[index]])
    )
  })
  patchwork::wrap_plots(panels, ncol = 1) +
    patchwork::plot_annotation(
      title = title,
      subtitle = subtitle,
      caption = paste(
        "One pooled 10x library per condition; cell-level statistics are descriptive.",
        "Attenuation is an effect-size criterion, not a formal interaction test."
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = 0),
        plot.subtitle = ggplot2::element_text(size = 9.3, color = "#333333", hjust = 0),
        plot.caption = ggplot2::element_text(size = 8, color = "#555555", hjust = 0)
      )
    )
}

figure13 <- make_figure(
  myeloid_levels,
  "Lenvatinib-associated transcriptional changes jointly attenuated across major myeloid lineages",
  "LEN DEG intersected with matched AMG44 and myeloid Eif2ak3-deletion attenuation; GO:BP uses each lineage-specific tested universe"
)
figure17 <- make_figure(
  lymphoid_levels,
  "Lenvatinib-associated transcriptional changes and downstream lymphoid ecosystem remodeling",
  "Joint attenuation in AMG44-treated and myeloid Eif2ak3-deficient backgrounds; lymphoid KO effects are not interpreted as cell intrinsic"
)

save_plot_atomic(figure13, output_paths$fig13_png, width = 16.6, height = 31.0)
save_plot_atomic(figure13, output_paths$fig13_pdf, width = 16.6, height = 31.0)
save_plot_atomic(figure17, output_paths$fig17_png, width = 16.6, height = 31.0)
save_plot_atomic(figure17, output_paths$fig17_pdf, width = 16.6, height = 31.0)

same_universe <- de_all |>
  dplyr::count(lineage, contrast, name = "n_tested") |>
  dplyr::group_by(lineage) |>
  dplyr::summarise(pass = dplyr::n_distinct(n_tested) == 1L, .groups = "drop") |>
  dplyr::pull(pass) |>
  all()
center_identity <- all(
  membership$joint_center ==
    (membership$LEN_DEG & membership$AMG_attenuation & membership$KO_attenuation)
)
attenuation_opposes_wt <- all(!effect_table$AMG_attenuated | effect_table$AMG_opposing_DiD) &&
  all(!effect_table$KO_attenuated | effect_table$KO_opposing_DiD)
center_genes_are_deg <- all(!membership$joint_center | membership$LEN_DEG)
go_display_is_simplified <- all(
  paste(go_display$lineage, go_display$direction, go_display$ID) %in%
    paste(go_simplified$lineage, go_simplified$direction, go_simplified$ID)
)
figure_sizes <- file.info(unlist(output_paths[c(
  "fig13_png", "fig13_pdf", "fig17_png", "fig17_pdf"
)]))$size

qa <- data.frame(
  check = c(
    "input_md5_contract", "canonical_cell_count", "canonical_feature_count",
    "six_fixed_groups", "eight_fixed_lineages", "all_lineage_groups_ge_30",
    "twenty_four_DE_comparisons", "same_universe_within_lineage",
    "explicit_LogNormalize_FC_formula", "BH_range", "finite_effects", "center_set_identity",
    "attenuation_opposes_WT", "center_genes_are_LEN_DEG",
    "sixteen_Venn_rows", "sixteen_GO_status_rows",
    "GO_only_with_at_least_10_mappable", "display_terms_are_significant",
    "display_terms_are_simplified", "display_term_cap_six",
    "four_final_figure_files_nonempty", "pooled_library_caveat"
  ),
  observed = c(
    input_md5,
    ncol(seu),
    nrow(counts),
    paste(group_levels, collapse = ","),
    paste(lineage_levels, collapse = ","),
    min(cell_counts$n_cells),
    nrow(de_all |> dplyr::distinct(lineage, contrast)),
    same_universe,
    max(fc_qa$max_abs_formula_error),
    all(is.na(de_all$p_adj_BH) | (de_all$p_adj_BH >= 0 & de_all$p_adj_BH <= 1)),
    all(is.finite(de_all$avg_log2FC)),
    center_identity,
    attenuation_opposes_wt,
    center_genes_are_deg,
    nrow(venn_counts),
    nrow(go_status),
    all(go_status$gate != "GO_ELIGIBLE" | go_status$n_mappable_center_genes >= 10L),
    all(is.finite(go_display$p.adjust) & go_display$p.adjust < 0.05),
    go_display_is_simplified,
    max(go_display |> dplyr::count(lineage, direction) |> dplyr::pull(n)),
    min(figure_sizes),
    TRUE
  ),
  expected = c(
    ifelse(is.na(config$expected_md5), "not enforced", config$expected_md5),
    ifelse(identical(input_md5, canonical_md5), 38774L, "not enforced"),
    ifelse(identical(input_md5, canonical_md5), 33696L, "not enforced"),
    paste(group_levels, collapse = ","),
    paste(lineage_levels, collapse = ","),
    ">=30", 24L, TRUE, "<=1e-12", TRUE, TRUE, TRUE, TRUE, TRUE,
    16L, 16L, TRUE, TRUE, TRUE, "<=6", ">1000 bytes", TRUE
  ),
  stringsAsFactors = FALSE
)
qa$status <- ifelse(c(
  is.na(config$expected_md5) || identical(tolower(input_md5), tolower(config$expected_md5)),
  !identical(input_md5, canonical_md5) || ncol(seu) == 38774L,
  !identical(input_md5, canonical_md5) || nrow(counts) == 33696L,
  setequal(unique(as.character(meta$group)), group_levels),
  setequal(unique(stats::na.omit(meta$lineage)), lineage_levels),
  all(cell_counts$n_cells >= min_cells_per_group),
  nrow(de_all |> dplyr::distinct(lineage, contrast)) == 24L,
  same_universe,
  all(is.finite(fc_qa$max_abs_formula_error)) &&
    max(fc_qa$max_abs_formula_error) <= 1e-12,
  all(is.na(de_all$p_adj_BH) | (de_all$p_adj_BH >= 0 & de_all$p_adj_BH <= 1)),
  all(is.finite(de_all$avg_log2FC)),
  center_identity,
  attenuation_opposes_wt,
  center_genes_are_deg,
  nrow(venn_counts) == 16L,
  nrow(go_status) == 16L,
  all(go_status$gate != "GO_ELIGIBLE" | go_status$n_mappable_center_genes >= 10L),
  all(is.finite(go_display$p.adjust) & go_display$p.adjust < 0.05),
  go_display_is_simplified,
  max(go_display |> dplyr::count(lineage, direction) |> dplyr::pull(n)) <= 6L,
  all(figure_sizes > 1000),
  TRUE
), "PASS", "FAIL")
if (any(qa$status == "FAIL")) {
  stop("Step7 QA failed: ", paste(qa$check[qa$status == "FAIL"], collapse = ", "), call. = FALSE)
}

analysis_contract <- data.frame(
  field = c(
    "script_version", "statistical_unit", "inference_scope", "groups",
    "lineages", "RNA_layer", "minimum_cells_per_lineage_group",
    "common_universe_rule", "DE_method", "DE_fold_change_mean_function",
    "DE_logfc_threshold_during_test",
    "DE_min_pct_during_test", "DEG_threshold", "multiple_testing",
    "AMG_DiD", "KO_DiD", "attenuation_rule", "joint_center",
    "GO_ontology", "GO_universe", "GO_adjustment", "GO_simplify",
    "GO_display_selection",
    "final_figures", "excluded_outputs"
  ),
  value = c(
    script_version,
    "one pooled 10x library per condition",
    paste(
      "Cell-level Wilcoxon/BH values are descriptive distributional summaries;",
      "not biological-replicate inference."
    ),
    paste(group_levels, collapse = ", "),
    paste(lineage_levels, collapse = ", "),
    "RNA log-normalized data; counts used only for detection universe",
    min_cells_per_group,
    "detected in >=10% of cells in any of six groups; same genes in WT/AMG/KO",
    "Seurat FindMarkers Wilcoxon",
    "explicit log2((rowSums(expm1(x))+1)/ncol(x)); object command provenance ignored",
    0,
    0,
    "absolute average log2FC >0.25 and within-comparison BH-FDR <0.05",
    "BH within each lineage x comparison",
    "(LENnAMG-CTRnAMG)-(LEN-CTR)",
    "(LENnKO-CTRnKO)-(LEN-CTR)",
    "opposing abs(DiD)>0.25 plus >=50% reduction or direction reversal",
    "LEN DEG AND AMG attenuation AND KO attenuation in matching LEN direction",
    "GO Biological Process",
    "lineage-specific tested gene universe",
    "BH; adjusted P <0.05; no threshold relaxation",
    "Wang semantic similarity, cutoff 0.70, select_fun=min",
    paste(
      "up to six simplified significant terms; prespecified terms have priority only",
      "when still significant, then fill by adjusted P, Count, and GO ID"
    ),
    "Figure 13 myeloid Venn+GO; Figure 17 lymphoid Venn+GO",
    "full-CD45 UMAP/composition; volcano; LEN-only GO; pDC/Baso DEG; other exploratory displays"
  ),
  stringsAsFactors = FALSE
)

input_manifest <- data.frame(
  role = "canonical_CD45_Seurat_RDS",
  path = normalizePath(config$input, winslash = "/", mustWork = TRUE),
  bytes = file.info(config$input)$size,
  md5 = input_md5,
  cells = ncol(seu),
  RNA_features = nrow(counts),
  stringsAsFactors = FALSE
)

write_csv_atomic(analysis_contract, output_paths$contract)
write_csv_atomic(input_manifest, output_paths$input_manifest)
write_csv_atomic(cell_counts, output_paths$cell_counts)
write_csv_atomic(gene_universe, output_paths$universe)
write_csv_atomic(de_all, output_paths$de_long)
write_csv_atomic(de_wide, output_paths$de_wide)
write_csv_atomic(fc_qa, output_paths$fc_qa)
write_csv_atomic(effect_table, output_paths$classification)
write_csv_atomic(membership, output_paths$membership)
write_csv_atomic(venn_counts, output_paths$venn_counts)
write_csv_atomic(center_genes, output_paths$center_genes)
write_csv_atomic(go_status, output_paths$go_status)
write_csv_atomic(go_mapping, output_paths$go_mapping)
write_csv_atomic(go_full, output_paths$go_full)
write_csv_atomic(go_simplified, output_paths$go_simplified)
write_csv_atomic(go_display, output_paths$go_display)
write_csv_atomic(qa, output_paths$qa)

summary_lines <- c(
  "# Step7 lineage-resolved CD45+ scRNA-seq analysis",
  "",
  paste0("- Script version: `", script_version, "`."),
  paste0("- Input MD5: `", input_md5, "`."),
  paste0("- Fixed groups: ", paste(group_levels, collapse = ", "), "."),
  paste0("- Fixed lineages: ", paste(lineage_levels, collapse = ", "), "."),
  "- Three matched-background effects were recalculated from the canonical Seurat object.",
  paste0(
    "- Fold change used an explicit LogNormalize mean function independent of stored Seurat commands; ",
    "maximum manual-formula error was ", format(max(fc_qa$max_abs_formula_error), scientific = TRUE), "."
  ),
  "- BH correction was recalculated within each lineage and comparison.",
  "- GO:BP used the corresponding lineage-specific tested universe and semantic simplify().",
  "- Only Figures 13 and 17 were generated.",
  "",
  "## Joint-center counts",
  "",
  paste0(
    "- ", venn_counts$lineage, " ", venn_counts$direction,
    ": ", venn_counts$joint_center_n, " / ", venn_counts$LEN_DEG_n,
    " LEN DEG jointly attenuated."
  ),
  "",
  "## Inference boundary",
  "",
  paste(
    "Each condition is represented by one pooled 10x library. Cell-level Wilcoxon",
    "tests and BH-adjusted P values are exploratory distributional summaries and",
    "do not provide treatment-level biological-replicate inference. Attenuation is",
    "a prespecified effect-size criterion, not a formal interaction or rescue test."
  )
)
write_lines_atomic(summary_lines, output_paths$report)
write_lines_atomic(capture.output(sessionInfo()), output_paths$session)

manifest_inputs <- unlist(output_paths[setdiff(names(output_paths), "output_manifest")])
manifest <- data.frame(
  role = names(manifest_inputs),
  path = normalizePath(manifest_inputs, winslash = "/", mustWork = TRUE),
  bytes = file.info(manifest_inputs)$size,
  md5 = unname(tools::md5sum(manifest_inputs)),
  stringsAsFactors = FALSE
)
write_csv_atomic(manifest, output_paths$output_manifest)

cat("STEP7_PASS\n")
cat("input_md5\t", input_md5, "\n", sep = "")
cat("output_dir\t", normalizePath(config$outdir, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("Neu_up_joint\t", venn_counts$joint_center_n[
  venn_counts$lineage == "Neu" & venn_counts$direction == "Up"
], "\n", sep = "")
cat("Neu_down_joint\t", venn_counts$joint_center_n[
  venn_counts$lineage == "Neu" & venn_counts$direction == "Down"
], "\n", sep = "")

