#!/usr/bin/env Rscript

# Reproduces the reported neutrophil query-native UMAP, five-cluster analysis,
# Ng label display, Xue tumour-state transfer, and Xue core-state maturation score.
# One pooled 10x library represents each condition; compositions and score
# distributions are descriptive and are not treatment-level statistical tests.

main <- function() {

# %% 00 - command-line parsing, packages, and immutable run contract

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1]]) || !nzchar(x[[1]])) y else x
}

parse_cli <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    token <- x[[i]]
    if (!startsWith(token, "--")) {
      stop("Unexpected positional argument: ", token, call. = FALSE)
    }
    token <- substring(token, 3L)
    if (grepl("=", token, fixed = TRUE)) {
      key <- sub("=.*$", "", token)
      value <- sub("^[^=]*=", "", token)
    } else {
      key <- token
      if (i == length(x) || startsWith(x[[i + 1L]], "--")) {
        value <- "true"
      } else {
        i <- i + 1L
        value <- x[[i]]
      }
    }
    key <- gsub("-", "_", key, fixed = TRUE)
    if (key %in% names(out)) stop("Duplicate argument: --", key, call. = FALSE)
    out[[key]] <- value
    i <- i + 1L
  }
  out
}

as_bool <- function(x, name) {
  value <- tolower(as.character(x))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("--", name, " must be true or false.", call. = FALSE)
}

resolve_project_path <- function(path, project_root) {
  is_absolute <- grepl("^([A-Za-z]:[/\\\\]|[/\\\\]{2}|/)", path)
  if (is_absolute) path else file.path(project_root, path)
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
prefix <- "Step6_Len_Neu_query_native_display"
project_root <- normalizePath(
  cli$project_root %||% ".", winslash = "/", mustWork = TRUE
)
query_file <- resolve_project_path(
  cli$query %||% file.path("results", "03_Len_Neu_NgMapping.rds"),
  project_root
)
xue_file <- resolve_project_path(
  cli$xue_reference %||%
    file.path("data", "reference", "seu_mouse_LiverCancer.anno.rds"),
  project_root
)
output_root <- resolve_project_path(
  cli$output_root %||% file.path("results", "reproducible_runs"),
  project_root
)
run_id <- cli$run_id %||% prefix
workers <- as.integer(cli$workers %||% "8")
hash_inputs <- as_bool(cli$hash_inputs %||% "true", "hash-inputs")

if (!is.finite(workers) || workers < 1L) {
  stop("--workers must be a positive integer.", call. = FALSE)
}
if (!grepl("^[A-Za-z0-9._-]+$", run_id)) {
  stop("--run-id may contain only letters, numbers, dot, underscore, and hyphen.")
}

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "tibble",
  "ggplot2", "patchwork", "ggridges", "scales", "presto", "UCell",
  "BiocParallel"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}
if (utils::packageVersion("Seurat") < "5.1.0") {
  stop("This script requires Seurat >= 5.1.0.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggridges)
})

options(stringsAsFactors = FALSE)

run_dir <- file.path(output_root, run_id)
if (file.exists(run_dir)) {
  stop("Refusing to overwrite existing run directory: ", run_dir, call. = FALSE)
}
dirs <- list(
  run = run_dir,
  figures = file.path(run_dir, "figures"),
  tables = file.path(run_dir, "tables"),
  objects = file.path(run_dir, "objects"),
  logs = file.path(run_dir, "logs")
)
for (path in dirs) dir.create(path, recursive = TRUE, showWarnings = FALSE)

run_complete <- FALSE
status_file <- file.path(dirs$logs, "RUN_STATUS.txt")
writeLines(
  c("status=RUNNING", paste0("started=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))),
  status_file
)
on.exit({
  if (!run_complete) {
    writeLines(
      c(
        "status=INCOMPLETE_OR_FAILED",
        paste0("updated=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        "Inspect the console/log and rerun with a new --run-id after correction."
      ),
      status_file
    )
  }
}, add = TRUE)

log_file <- file.path(dirs$logs, paste0(prefix, "__console.log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)
log_closed <- FALSE
close_log <- function(silent = FALSE) {
  if (log_closed) return(invisible(NULL))
  close_connections <- function() {
    if (sink.number(type = "message") != 2L) sink(type = "message")
    if (sink.number(type = "output") > 0L) sink(type = "output")
    if (isOpen(log_con)) close(log_con)
  }
  if (silent) {
    try(close_connections(), silent = TRUE)
  } else {
    close_connections()
  }
  log_closed <<- TRUE
  invisible(NULL)
}
on.exit(close_log(silent = TRUE), add = TRUE)

stage <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste0(..., collapse = ""))
}

write_csv <- function(x, stem) {
  path <- file.path(dirs$tables, paste0(prefix, "__", stem, ".csv"))
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  path
}

save_rds <- function(x, stem, compress = TRUE) {
  path <- file.path(dirs$objects, paste0(prefix, "__", stem, ".rds"))
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  saveRDS(x, path, compress = compress)
  path
}

save_plot_both <- function(plot, stem, width, height) {
  base <- file.path(dirs$figures, paste0(prefix, "__", stem))
  png_file <- paste0(base, ".png")
  pdf_file <- paste0(base, ".pdf")
  if (file.exists(png_file) || file.exists(pdf_file)) {
    stop("Refusing to overwrite plot output: ", stem)
  }
  ggsave(png_file, plot = plot, width = width, height = height,
         units = "in", dpi = 320, bg = "white")
  ggsave(pdf_file, plot = plot, width = width, height = height,
         units = "in", bg = "white")
  c(png_file, pdf_file)
}

rank01 <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(out)
  if (sum(ok) == 1L || length(unique(x[ok])) == 1L) {
    out[ok] <- 0.5
  } else {
    out[ok] <- (rank(x[ok], ties.method = "average") - 1) / (sum(ok) - 1)
  }
  out
}

rescale01 <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(out)
  observed <- range(x[ok])
  if (diff(observed) == 0) {
    out[ok] <- 0.5
  } else {
    out[ok] <- (x[ok] - observed[[1]]) / diff(observed)
  }
  out
}

technical_hvg_blacklist <- function(genes) {
  grepl(
    paste0(
      "^(mt-|Rpl|Rps|Hba|Hbb|Gm[0-9]+|Xist$|Ddx3y$|Eif2s3y$|",
      "Kdm5d$|Uty$|Malat1$|Pf4$|Ppbp$)"
    ),
    genes,
    ignore.case = TRUE
  )
}

select_balanced_hvgs <- function(counts, groups, group_levels,
                                 nfeatures_per_group = 4000L,
                                 target_n = 3000L) {
  within_group <- vector("list", length(group_levels))
  names(within_group) <- group_levels
  for (group_name in group_levels) {
    cells <- colnames(counts)[groups[colnames(counts)] == group_name]
    if (length(cells) < 50L) stop("Too few cells in group ", group_name)
    stage("Selecting within-library HVGs: ", group_name, " (n=", length(cells), ")")
    temp <- CreateSeuratObject(
      counts = counts[, cells, drop = FALSE], min.cells = 0, min.features = 0
    )
    temp <- NormalizeData(temp, normalization.method = "LogNormalize",
                          scale.factor = 10000, verbose = FALSE)
    temp <- FindVariableFeatures(
      temp, selection.method = "vst",
      nfeatures = min(nfeatures_per_group, nrow(temp)), verbose = FALSE
    )
    genes <- VariableFeatures(temp)
    within_group[[group_name]] <- tibble(
      group = group_name,
      gene = genes,
      within_group_rank = seq_along(genes)
    )
    rm(temp)
    gc()
  }
  within_long <- bind_rows(within_group)
  ranking <- within_long |>
    filter(!technical_hvg_blacklist(gene)) |>
    group_by(gene) |>
    summarise(
      n_libraries = n_distinct(group),
      median_rank = median(within_group_rank),
      mean_rank = mean(within_group_rank),
      best_rank = min(within_group_rank),
      .groups = "drop"
    ) |>
    arrange(desc(n_libraries), median_rank, mean_rank, best_rank, gene) |>
    mutate(selected = row_number() <= target_n)
  selected <- ranking$gene[ranking$selected]
  if (length(selected) != target_n) {
    stop("Balanced HVG selection returned ", length(selected), " rather than ", target_n)
  }
  list(selected = selected, ranking = ranking, within_group = within_long)
}

lognorm_mean_fxn <- function(x) {
  log2((Matrix::rowSums(expm1(x)) + 1) / ncol(x))
}

join_assay_layers <- function(object, assay) {
  count_layers <- Layers(object[[assay]], search = "^counts")
  data_layers <- Layers(object[[assay]], search = "^data")
  if (length(count_layers) > 1L || length(data_layers) > 1L) {
    object[[assay]] <- JoinLayers(object[[assay]])
  }
  object
}

group_levels <- c("CTR", "LEN", "CTRnAMG", "LENnAMG", "CTRnKO", "LENnKO")
ng_stage_levels <- c(
  "preNeu", "IMM 1", "IMM 2", "MAT 1", "MAT 2", "MAT 3", "MAT 5",
  "T1", "T2", "T3"
)
xue_tumour_state_levels <- c(
  "mNeu_07_Actg1", "mNeu_08_Ltf", "mNeu_09_Apoa2",
  "mNeu_10_Ifit1", "mNeu_11_Spp1", "mNeu_12_Ccl4"
)
xue_core_early_states <- c("mNeu_01_Ngp", "mNeu_02_Mmp8")
xue_core_mature_states <- c(
  "mNeu_03_Pabpc1", "mNeu_05_Gm2a", "mNeu_06_Marco", "mNeu_07_Actg1"
)
xue_core_state_levels <- c(xue_core_early_states, xue_core_mature_states)

query_umap_seed <- 123L
# Fixed seeds reproduce the reported reference mapping and cluster labels.
mapping_seed <- 20260813L
cluster_seed <- 20260815L
set.seed(query_umap_seed)

input_paths <- c(query = query_file, xue_reference = xue_file)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop("Missing input(s):\n", paste(missing_inputs, collapse = "\n"), call. = FALSE)
}
input_info <- file.info(input_paths)
if (any(is.na(input_info$size)) || any(input_info$size <= 0)) {
  stop("At least one input file is empty or unreadable.")
}
input_manifest <- tibble(
  input = names(input_paths),
  path = unname(input_paths),
  size_bytes = input_info$size,
  modified = format(input_info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
  md5 = if (hash_inputs) unname(tools::md5sum(input_paths)) else NA_character_
)
input_manifest_file <- write_csv(input_manifest, "input_manifest")
print(input_manifest)

# %% 01 - load and validate the upstream Ng-mapped study Neu query

stage("Loading upstream Ng-mapped study neutrophil query.")
query <- readRDS(query_file)
if (!inherits(query, "Seurat")) stop("Query is not a Seurat object.")
required_query_meta <- c("group", "CellType_l1", "Ng_stage")
if (!all(required_query_meta %in% colnames(query[[]]))) {
  stop("Query lacks metadata: ",
       paste(setdiff(required_query_meta, colnames(query[[]])), collapse = ", "))
}
if (!all(c("RNA", "integrated") %in% Assays(query))) {
  stop("Query must contain RNA and integrated assays.")
}
if (ncol(query) != 14111L || !all(as.character(query$CellType_l1) == "Neu")) {
  stop("Expected the frozen 14,111-cell Neu query from upstream Step3.")
}
if (!setequal(unique(as.character(query$group)), group_levels)) {
  stop("Query does not contain the expected six pooled libraries.")
}
if (!setequal(unique(as.character(query$Ng_stage)), ng_stage_levels)) {
  stop("Unexpected Ng_stage levels in the upstream mapped query.")
}
query$group <- factor(as.character(query$group), levels = group_levels)
query$Ng_stage <- factor(as.character(query$Ng_stage), levels = ng_stage_levels)
query <- join_assay_layers(query, "RNA")
DefaultAssay(query) <- "RNA"
if (!"counts" %in% Layers(query[["RNA"]])) stop("Query RNA counts layer is absent.")
query <- NormalizeData(
  query, assay = "RNA", normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = FALSE
)

# %% 02 - independently calculate the query-native Neu UMAP

stage("Building the independently calculated query-native Neu UMAP.")
query_counts <- LayerData(query, assay = "RNA", layer = "counts")
balanced_hvg <- select_balanced_hvgs(
  counts = query_counts,
  groups = setNames(as.character(query$group), colnames(query)),
  group_levels = group_levels,
  nfeatures_per_group = 4000L,
  target_n = 3000L
)
query <- ScaleData(
  query, assay = "RNA", features = balanced_hvg$selected, verbose = FALSE
)
query <- RunPCA(
  query, assay = "RNA", features = balanced_hvg$selected, npcs = 50,
  reduction.name = "query.neu.pca", reduction.key = "QUERYNEUPC_",
  seed.use = query_umap_seed, verbose = FALSE
)
query <- RunUMAP(
  query, reduction = "query.neu.pca", dims = 1:30,
  n.neighbors = 50, min.dist = 0.35, metric = "cosine",
  reduction.name = "query.neu.umap", reduction.key = "QUERYNEUUMAP_",
  seed.use = query_umap_seed, verbose = TRUE
)
query_native_umap_before_mapping <- Embeddings(query, "query.neu.umap")
stopifnot(
  identical(rownames(query_native_umap_before_mapping), colnames(query)),
  all(is.finite(query_native_umap_before_mapping))
)
rm(query_counts)
gc()

# %% 03 - load Xue reference; construct the retained core-state marker axis

stage("Loading the canonical Xue/Zhang mouse liver-tumour reference.")
xue_all <- readRDS(xue_file)
if (!inherits(xue_all, "Seurat") || !"RNA" %in% Assays(xue_all)) {
  stop("Xue input is not the expected Seurat object with an RNA assay.")
}
required_xue_meta <- c("Sample", "celltype", "clusters", "tissue")
if (!all(required_xue_meta %in% colnames(xue_all[[]]))) {
  stop("Xue reference lacks required metadata: ",
       paste(setdiff(required_xue_meta, colnames(xue_all[[]])), collapse = ", "))
}
xue_mneu_all <- subset(xue_all, subset = celltype == "mNeu")
if (ncol(xue_mneu_all) != 17780L) {
  stop("Expected 17,780 mouse mNeu cells in the canonical Xue object.")
}
xue_all_mneu_n <- ncol(xue_mneu_all)
rm(xue_all)
gc()
xue_mneu_all <- join_assay_layers(xue_mneu_all, "RNA")
DefaultAssay(xue_mneu_all) <- "RNA"
xue_mneu_all <- NormalizeData(
  xue_mneu_all, assay = "RNA", normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = FALSE
)
xue_all_state_levels <- sort(unique(as.character(xue_mneu_all$clusters)))
if (!all(c(xue_tumour_state_levels, xue_core_state_levels) %in% xue_all_state_levels)) {
  stop("Expected Xue mNeu states are absent from the canonical reference.")
}

stage("Deriving Xue core-state marker sets from all 17,780 mouse mNeu cells.")
xue_core_marker_stats <- presto::wilcoxauc(
  xue_mneu_all, group_by = "clusters", assay = "data", seurat_assay = "RNA"
) |>
  as_tibble() |>
  mutate(group = as.character(group))
xue_core_marker_reference <- xue_core_marker_stats |>
  filter(
    group %in% xue_core_state_levels,
    padj < 0.05,
    logFC >= 0.25,
    pct_in >= 0.10,
    auc >= 0.55
  ) |>
  group_by(group) |>
  arrange(desc(auc), desc(logFC), .by_group = TRUE) |>
  slice_head(n = 30L) |>
  ungroup()
xue_core_reference_counts <- xue_core_marker_reference |>
  count(group, name = "n_reference_markers")
if (!setequal(xue_core_reference_counts$group, xue_core_state_levels) ||
    any(xue_core_reference_counts$n_reference_markers != 30L)) {
  stop("Failed to derive exactly 30 Xue core markers for each retained state.")
}
rm(xue_core_marker_stats)
gc()

# %% 04 - select the tumour-only pTMC/pTMK Xue states 07-12

xue_mneu_meta <- xue_mneu_all[[]] |>
  rownames_to_column("cell") |>
  mutate(
    Sample = as.character(Sample),
    tissue = as.character(tissue),
    clusters = as.character(clusters),
    is_pTMC_or_pTMK = grepl("pTMC|pTMK", Sample, ignore.case = TRUE)
  )
xue_reference_audit <- xue_mneu_meta |>
  count(Sample, tissue, clusters, is_pTMC_or_pTMK, name = "n_cells")

target_mask <- xue_mneu_meta$tissue == "Tumor" &
  xue_mneu_meta$clusters %in% xue_tumour_state_levels
if (any(target_mask & !xue_mneu_meta$is_pTMC_or_pTMK)) {
  stop("Tumour-state 07-12 cells include samples outside pTMC/pTMK.")
}
xue_tumour_cells <- xue_mneu_meta$cell[target_mask & xue_mneu_meta$is_pTMC_or_pTMK]
xue_ref <- subset(xue_mneu_all, cells = xue_tumour_cells)
if (ncol(xue_ref) != 8297L) {
  stop("Expected 8,297 pTMC/pTMK tumour-derived Xue mNeu cells; observed ", ncol(xue_ref))
}
xue_ref$Xue_state <- factor(
  as.character(xue_ref$clusters), levels = xue_tumour_state_levels
)
if (anyNA(xue_ref$Xue_state) || !all(as.character(xue_ref$tissue) == "Tumor")) {
  stop("Invalid Xue tumour-reference labels.")
}
Idents(xue_ref) <- "Xue_state"
rm(xue_mneu_all)
gc()

# %% 05 - rebuild Xue reference PCA/UMAP model and transfer states

stage("Rebuilding the Xue tumour-mNeu mapping reference.")
DefaultAssay(xue_ref) <- "RNA"
xue_ref <- NormalizeData(
  xue_ref, assay = "RNA", normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = FALSE
)
xue_ref <- FindVariableFeatures(
  xue_ref, assay = "RNA", selection.method = "vst", nfeatures = 2000,
  verbose = FALSE
)
mapping_features <- intersect(VariableFeatures(xue_ref), rownames(query))
if (length(mapping_features) < 1500L) {
  stop("Fewer than 1,500 Xue reference HVGs are shared with the query.")
}
xue_ref <- ScaleData(
  xue_ref, assay = "RNA", features = mapping_features, verbose = FALSE
)
xue_ref <- RunPCA(
  xue_ref, assay = "RNA", features = mapping_features, npcs = 50,
  reduction.name = "xue.pca", reduction.key = "XUEPC_",
  seed.use = mapping_seed, verbose = FALSE
)
xue_ref <- RunUMAP(
  xue_ref, reduction = "xue.pca", dims = 1:30,
  n.neighbors = 30, min.dist = 0.30, metric = "cosine",
  reduction.name = "xue.umap", reduction.key = "XUEUMAP_",
  return.model = TRUE, seed.use = mapping_seed, verbose = TRUE
)

stage("Finding Xue-to-query transfer anchors and running MapQuery.")
set.seed(mapping_seed)
xue_anchors <- FindTransferAnchors(
  reference = xue_ref,
  query = query,
  normalization.method = "LogNormalize",
  reference.assay = "RNA",
  query.assay = "RNA",
  reduction = "pcaproject",
  reference.reduction = "xue.pca",
  features = mapping_features,
  dims = 1:30,
  verbose = TRUE
)
if (nrow(xue_anchors@anchors) == 0L) stop("No Xue-query anchors were found.")

query <- MapQuery(
  anchorset = xue_anchors,
  query = query,
  reference = xue_ref,
  refdata = list(Xue_state = "Xue_state"),
  new.reduction.name = "xue.ref.pca",
  reference.reduction = "xue.pca",
  reduction.model = "xue.umap",
  store.weights = FALSE,
  projectumap.args = list(
    reduction.name = "xue.ref.umap",
    reduction.key = "XUEREFUMAP_"
  ),
  verbose = TRUE
)
query_native_umap_after_mapping <- Embeddings(query, "query.neu.umap")
native_umap_max_abs_difference <- max(abs(
  query_native_umap_after_mapping - query_native_umap_before_mapping
))
if (native_umap_max_abs_difference > 1e-12) {
  stop("MapQuery altered the independently calculated query-native UMAP.")
}

if (!all(c("predicted.Xue_state", "Ng_stage") %in% colnames(query[[]])) ||
    !"prediction.score.Xue_state" %in% Assays(query)) {
  stop("MapQuery did not return the expected Xue labels and score assay.")
}
query$Xue_state <- factor(
  as.character(query$predicted.Xue_state), levels = xue_tumour_state_levels
)
if (anyNA(query$Xue_state)) stop("Xue label transfer returned missing labels.")

xue_score_matrix <- as.matrix(LayerData(
  query, assay = "prediction.score.Xue_state", layer = "data"
))
xue_score_names <- gsub("-", "_", rownames(xue_score_matrix), fixed = TRUE)
if (setequal(xue_score_names, xue_tumour_state_levels)) {
  rownames(xue_score_matrix) <- xue_score_names
}
if (!setequal(rownames(xue_score_matrix), xue_tumour_state_levels) ||
    !setequal(colnames(xue_score_matrix), colnames(query))) {
  stop("Xue transferred score matrix is not aligned to the query.")
}
xue_score_matrix <- xue_score_matrix[
  xue_tumour_state_levels, colnames(query), drop = FALSE
]
score_order <- apply(xue_score_matrix, 2, sort, decreasing = TRUE)
prediction_from_score <- rownames(xue_score_matrix)[
  max.col(t(xue_score_matrix), ties.method = "first")
]
query$Xue_prediction_max <- score_order[1, ]
query$Xue_prediction_second <- score_order[2, ]
query$Xue_prediction_margin <- score_order[1, ] - score_order[2, ]
label_score_agreement <- mean(
  prediction_from_score == as.character(query$Xue_state)
)
if (label_score_agreement < 0.999 ||
    any(query$Xue_prediction_max < 0 | query$Xue_prediction_max > 1 + 1e-8) ||
    any(query$Xue_prediction_margin < -1e-8)) {
  stop("Transferred Xue labels or confidence scores failed validation.")
}

# %% 06 - calculate only the retained Xue core-state maturation score

stage("Calculating the retained Xue core-state maturation score with UCell.")
xue_core_marker_manifest <- xue_core_marker_reference |>
  mutate(
    mapped_feature = rownames(query)[match(toupper(feature), toupper(rownames(query)))],
    present_in_query = !is.na(mapped_feature),
    core_role = if_else(
      group %in% xue_core_early_states, "early_01_02", "mature_like_03_05_06_07"
    )
  )
xue_core_marker_used <- xue_core_marker_manifest |>
  filter(present_in_query) |>
  distinct(group, mapped_feature, .keep_all = TRUE)
xue_core_coverage <- xue_core_marker_used |>
  count(group, core_role, name = "n_genes_detected")
if (!setequal(xue_core_coverage$group, xue_core_state_levels) ||
    any(xue_core_coverage$n_genes_detected < 5L)) {
  stop("At least one Xue core state has fewer than five query-detected markers.")
}
xue_core_signatures <- split(
  xue_core_marker_used$mapped_feature, xue_core_marker_used$group
)
xue_core_signatures <- xue_core_signatures[xue_core_state_levels]
names(xue_core_signatures) <- paste0("XueCore_", xue_core_state_levels)

detected_cores <- parallel::detectCores(logical = FALSE)
if (!is.finite(detected_cores)) detected_cores <- 1L
score_workers <- max(1L, min(workers, as.integer(detected_cores), 8L))
ucell_bpparam <- if (.Platform$OS.type == "unix" && score_workers > 1L) {
  BiocParallel::MulticoreParam(
    workers = score_workers, progressbar = TRUE, stop.on.error = TRUE
  )
} else {
  BiocParallel::SerialParam(progressbar = TRUE, stop.on.error = TRUE)
}
query <- UCell::AddModuleScore_UCell(
  obj = query,
  features = xue_core_signatures,
  assay = "RNA",
  BPPARAM = ucell_bpparam,
  ncores = 1L,
  name = "_UCell"
)
xue_core_score_columns <- paste0(
  "XueCore_", xue_core_state_levels, "_UCell"
)
if (!all(xue_core_score_columns %in% colnames(query[[]]))) {
  stop("UCell did not return all Xue core-state score columns.")
}
xue_core_rank_matrix <- sapply(
  xue_core_score_columns,
  function(column) rank01(query[[column, drop = TRUE]])
)
colnames(xue_core_rank_matrix) <- xue_core_score_columns
early_columns <- paste0("XueCore_", xue_core_early_states, "_UCell")
mature_columns <- paste0("XueCore_", xue_core_mature_states, "_UCell")
xue_core_axis_raw <-
  rowMeans(xue_core_rank_matrix[, mature_columns, drop = FALSE]) -
  rowMeans(xue_core_rank_matrix[, early_columns, drop = FALSE])
query$Xue_core_maturation <- rescale01(xue_core_axis_raw)
if (anyNA(query$Xue_core_maturation) ||
    any(query$Xue_core_maturation < 0 | query$Xue_core_maturation > 1)) {
  stop("Xue core-state maturation score failed range validation.")
}

# %% 07 - de novo integrated-assay clustering at the retained resolution 0.3

stage("Running the retained de novo five-cluster Neu analysis.")
DefaultAssay(query) <- "integrated"
integrated_features <- VariableFeatures(query[["integrated"]])
if (length(integrated_features) < 500L) {
  stop("Integrated assay has fewer than 500 variable features.")
}
integrated_scaled <- LayerData(query, assay = "integrated", layer = "scale.data")
if (nrow(integrated_scaled) == 0L) {
  query <- ScaleData(
    query, assay = "integrated", features = integrated_features, verbose = FALSE
  )
}
rm(integrated_scaled)
gc()
query <- RunPCA(
  query, assay = "integrated", features = integrated_features, npcs = 30,
  reduction.name = "neu.integrated.pca",
  reduction.key = "NEUINTPC_", seed.use = cluster_seed, verbose = FALSE
)
query <- FindNeighbors(
  query, reduction = "neu.integrated.pca", dims = 1:20,
  k.param = 50, graph.name = c("neu_nn", "neu_snn"),
  verbose = TRUE
)
query <- FindClusters(
  query, graph.name = "neu_snn", resolution = 0.3,
  algorithm = 1, random.seed = cluster_seed, verbose = TRUE
)

raw_cluster <- as.character(Idents(query))
raw_levels <- sort(unique(as.integer(raw_cluster)))
if (length(raw_levels) != 5L || anyNA(raw_levels)) {
  stop(
    "Expected five de novo clusters at resolution 0.3; observed: ",
    paste(sort(unique(raw_cluster)), collapse = ", ")
  )
}
initial_lookup <- setNames(paste0("C", seq_along(raw_levels)), raw_levels)
initial_labels <- unname(initial_lookup[raw_cluster])
initial_levels <- unname(initial_lookup[as.character(raw_levels)])
query$Neu_cluster_res0_3_initial <- factor(
  initial_labels, levels = initial_levels
)

# This retained relabel changes neither clustering nor UMAP coordinates.
cluster_relabel <- tibble(
  old_cluster = c("C5", "C3", "C4", "C1", "C2"),
  new_cluster = paste0("C", 1:5),
  display_cluster = c(
    "C1_Ltf", "C2_Cxcr2", "C3_Ifit1", "C4_Olr1", "C5_Ccl4"
  ),
  requested_order = 1:5
)
if (!setequal(cluster_relabel$old_cluster, initial_levels)) {
  stop("Initial cluster identities do not match the frozen relabel contract.")
}
new_lookup <- setNames(cluster_relabel$new_cluster, cluster_relabel$old_cluster)
display_lookup <- setNames(
  cluster_relabel$display_cluster, cluster_relabel$new_cluster
)
numeric_labels <- unname(new_lookup[initial_labels])
display_labels <- unname(display_lookup[numeric_labels])
cluster_numeric_levels <- cluster_relabel$new_cluster
cluster_display_levels <- cluster_relabel$display_cluster
query$Neu_cluster_res0_3 <- factor(
  numeric_labels, levels = cluster_numeric_levels
)
query$Neu_cluster_res0_3_display <- factor(
  display_labels, levels = cluster_display_levels
)
Idents(query) <- "Neu_cluster_res0_3_display"

cluster_relabel_score_audit <- tibble(
  old_cluster = initial_labels,
  Xue_core_maturation = query$Xue_core_maturation
) |>
  group_by(old_cluster) |>
  summarise(
    n_cells = n(),
    Xue_core_mean_before_relabel = mean(Xue_core_maturation),
    Xue_core_median_before_relabel = median(Xue_core_maturation),
    .groups = "drop"
  )
cluster_relabel <- cluster_relabel |>
  left_join(cluster_relabel_score_audit, by = "old_cluster") |>
  arrange(requested_order) |>
  mutate(
    relabel_basis = paste0(
      "Retained Xue-core-informed order; graph clusters and query-native UMAP unchanged"
    )
  )

query_native_umap_after_clustering <- Embeddings(query, "query.neu.umap")
native_umap_cluster_max_abs_difference <- max(abs(
  query_native_umap_after_clustering - query_native_umap_before_mapping
))
if (native_umap_cluster_max_abs_difference > 1e-12) {
  stop("De novo clustering altered the query-native UMAP coordinates.")
}

# %% 08 - marker ranking and retained descriptive tables

stage("Ranking de novo cluster markers for annotation only.")
DefaultAssay(query) <- "RNA"
query <- join_assay_layers(query, "RNA")
Idents(query) <- "Neu_cluster_res0_3_display"
cluster_markers <- FindAllMarkers(
  query,
  assay = "RNA",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25,
  return.thresh = 0.05,
  mean.fxn = lognorm_mean_fxn,
  fc.name = "avg_log2FC",
  verbose = TRUE
) |>
  as_tibble() |>
  mutate(
    cluster = factor(as.character(cluster), levels = cluster_display_levels),
    analysis_scope = paste0(
      "cluster-marker ranking for annotation; not a treatment-level statistical test"
    )
  ) |>
  arrange(cluster, p_val_adj, desc(avg_log2FC), gene)
if (nrow(cluster_markers) == 0L ||
    !setequal(as.character(unique(cluster_markers$cluster)), cluster_display_levels)) {
  stop("FindAllMarkers did not return markers for all five retained clusters.")
}
cluster_markers_top20 <- cluster_markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 20L, with_ties = FALSE) |>
  arrange(cluster, desc(avg_log2FC), p_val_adj, gene) |>
  ungroup()

marker_categories <- list(
  `Immature / granules` = c(
    "S100a8", "S100a9", "Retnlg", "Ngp", "Camp", "Ltf", "Mmp8", "Lcn2"
  ),
  `Maturation / trafficking` = c("Csf3r", "Cxcr2", "Fpr1", "Fpr2", "Itgam"),
  Inflammatory = c("Il1b", "Cxcl2", "Ccl3", "Ccl4", "Spp1"),
  `IFN / ISG` = c("Ifit1", "Ifit3", "Isg15", "Stat1", "Cxcl10"),
  Regulatory = c("Cd274", "Arg1", "Olr1", "Il1rn", "Marco"),
  `ISR / TF` = c("Atf4", "Ddit3", "Ppp1r15a", "Cebpd")
)
marker_manifest <- bind_rows(lapply(names(marker_categories), function(category_name) {
  genes <- marker_categories[[category_name]]
  tibble(
    category = category_name,
    gene = genes,
    detected_in_query = genes %in% rownames(query)
  )
}))
if (nrow(marker_manifest) != sum(lengths(marker_categories))) {
  stop("Marker manifest row count is incomplete.")
}
marker_categories_detected <- lapply(
  marker_categories, intersect, y = rownames(query)
)
if (any(lengths(marker_categories_detected) < 3L)) {
  stop("A retained marker category has fewer than three detected genes.")
}

native_coordinates <- as.data.frame(Embeddings(query, "query.neu.umap")) |>
  rownames_to_column("cell")
colnames(native_coordinates)[2:3] <- c("UMAP_1", "UMAP_2")
per_cell <- query[[]] |>
  rownames_to_column("cell") |>
  transmute(
    cell,
    group = factor(as.character(group), levels = group_levels),
    Neu_cluster_initial = factor(
      as.character(Neu_cluster_res0_3_initial), levels = initial_levels
    ),
    Neu_cluster = factor(
      as.character(Neu_cluster_res0_3_display), levels = cluster_display_levels
    ),
    Ng_stage = factor(as.character(Ng_stage), levels = ng_stage_levels),
    Xue_state = factor(as.character(Xue_state), levels = xue_tumour_state_levels),
    Xue_prediction_max,
    Xue_prediction_second,
    Xue_prediction_margin,
    Xue_core_maturation
  ) |>
  left_join(native_coordinates, by = "cell")
if (nrow(per_cell) != 14111L || anyNA(per_cell)) {
  stop("Per-cell audit table is incomplete.")
}

cluster_sizes <- per_cell |>
  count(Neu_cluster, name = "n_cells", .drop = FALSE)
cluster_composition <- per_cell |>
  count(group, Neu_cluster, name = "n_cells", .drop = FALSE) |>
  group_by(group) |>
  mutate(percent_within_captured_Neu = 100 * n_cells / sum(n_cells)) |>
  ungroup()
maturation_summary <- per_cell |>
  group_by(Neu_cluster) |>
  summarise(
    n_cells = n(),
    mean = mean(Xue_core_maturation),
    median = median(Xue_core_maturation),
    q25 = quantile(Xue_core_maturation, 0.25),
    q75 = quantile(Xue_core_maturation, 0.75),
    score_definition = paste0(
      "Xue all-mNeu core markers: rank-mean(states 03/05/06/07) - ",
      "rank-mean(states 01/02), rescaled to 0-1"
    ),
    statistical_unit = paste0(
      "cells are descriptive observations from pooled 10x libraries; no treatment test"
    ),
    .groups = "drop"
  )
mapping_qc_by_group <- per_cell |>
  group_by(group) |>
  summarise(
    n_cells = n(),
    median_max_score = median(Xue_prediction_max),
    q25_max_score = quantile(Xue_prediction_max, 0.25),
    q75_max_score = quantile(Xue_prediction_max, 0.75),
    median_margin = median(Xue_prediction_margin),
    q25_margin = quantile(Xue_prediction_margin, 0.25),
    q75_margin = quantile(Xue_prediction_margin, 0.75),
    .groups = "drop"
  )
xue_composition <- per_cell |>
  count(group, Xue_state, name = "n_cells", .drop = FALSE) |>
  group_by(group) |>
  mutate(percent_within_captured_Neu = 100 * n_cells / sum(n_cells)) |>
  ungroup()
ng_xue_crosswalk <- per_cell |>
  count(Ng_stage, Xue_state, name = "n_cells", .drop = FALSE) |>
  group_by(Ng_stage) |>
  mutate(percent_within_Ng_stage = 100 * n_cells / sum(n_cells)) |>
  ungroup()
xue_reference_counts <- xue_ref[[]] |>
  rownames_to_column("cell") |>
  count(Sample, Xue_state, tissue, name = "n_cells", .drop = FALSE)
xue_numeric_crosswalk <- tibble(
  numeric_state = 7:12,
  object_label = xue_tumour_state_levels,
  paper_fig4_label = c(
    "mNeu_07_Actg1", "mNeu_08_Mmp8", "mNeu_09_Apoa2",
    "mNeu_10_Ifit1", "mNeu_11_Spp1", "mNeu_12_Ccl4"
  ),
  note = c(
    "object and paper suffix agree",
    "downloaded object/MOESM6 uses Ltf; paper Fig. 4/MOESM5 uses Mmp8",
    rep("object and paper suffix agree", 4)
  )
)
mapping_qc_overall <- tibble(
  metric = c(
    "query_cells", "Xue_all_mNeu_cells_for_core_markers",
    "Xue_tumour_reference_cells", "Xue_tumour_states",
    "shared_mapping_features", "mapping_dimensions", "anchor_pairs",
    "label_score_agreement", "median_max_score", "median_margin",
    "query_native_UMAP_change_after_mapping",
    "query_native_UMAP_change_after_clustering"
  ),
  value = as.character(c(
    ncol(query), xue_all_mneu_n, ncol(xue_ref),
    length(xue_tumour_state_levels), length(mapping_features), 30L,
    nrow(xue_anchors@anchors), label_score_agreement,
    median(query$Xue_prediction_max), median(query$Xue_prediction_margin),
    native_umap_max_abs_difference, native_umap_cluster_max_abs_difference
  )),
  interpretation = c(
    "frozen study Neu query", "source population used for marker discovery",
    "pTMC/pTMK tumour cells in states 07-12", "transferred categorical states",
    "Xue HVGs shared with query", "PCs 1-30", "Seurat transfer anchors",
    "argmax transferred score vs predicted label", "mapping QC", "mapping QC",
    "must equal zero", "must equal zero"
  )
)

# %% 09 - reported figure components

stage("Building the reported figure components.")
cluster_palette <- c(
  "C1_Ltf" = "#4477AA", "C2_Cxcr2" = "#EE6677",
  "C3_Ifit1" = "#228833", "C4_Olr1" = "#CCBB44",
  "C5_Ccl4" = "#66CCEE"
)
ng_palette <- c(
  "preNeu" = "#0072B2", "IMM 1" = "#56B4E9", "IMM 2" = "#88CCEE",
  "MAT 1" = "#44AA99", "MAT 2" = "#32806E", "MAT 3" = "#888888",
  "MAT 5" = "#B54C86", "T1" = "#D55E00", "T2" = "#F0CA42",
  "T3" = "#F79000"
)
xue_palette <- c(
  "mNeu_07_Actg1" = "#F4BA6F", "mNeu_08_Ltf" = "#EE7B1B",
  "mNeu_09_Apoa2" = "#CAB2D6", "mNeu_10_Ifit1" = "#6A3D9A",
  "mNeu_11_Spp1" = "#FFD700", "mNeu_12_Ccl4" = "#B15928"
)

umap_theme <- theme_classic(base_size = 10) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )
umap_x_range <- range(per_cell$UMAP_1)
umap_y_range <- range(per_cell$UMAP_2)
umap_x_range <- umap_x_range + c(-1, 1) * diff(umap_x_range) * 0.015
umap_y_range <- umap_y_range + c(-1, 1) * diff(umap_y_range) * 0.015
native_coord <- function() {
  coord_fixed(1, xlim = umap_x_range, ylim = umap_y_range, expand = FALSE)
}

p_marker_dotplot <- DotPlot(
  query,
  assay = "RNA",
  group.by = "Neu_cluster_res0_3_display",
  features = marker_categories_detected,
  dot.scale = 6,
  scale = TRUE
) +
  scale_color_gradientn(
    colours = c("white", "#FDE3D8", "#ED684E", "#D21E20", "#981917"),
    name = "Scaled average\nexpression"
  ) +
  labs(
    title = "A  Selected marker genes across de novo TAN clusters",
    subtitle = "Dot size: percentage expressed; color: scaled average expression"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.title = element_text(face = "bold")
  )
marker_dotplot_data <- as_tibble(p_marker_dotplot$data)

ctr_len <- per_cell |>
  filter(group %in% c("CTR", "LEN")) |>
  mutate(group = factor(as.character(group), levels = c("CTR", "LEN")))
p_cluster_ctr_len <- ggplot(
  ctr_len, aes(UMAP_1, UMAP_2, color = Neu_cluster)
) +
  geom_point(size = 0.42, alpha = 0.88, stroke = 0) +
  facet_wrap(~group, nrow = 1, labeller = as_labeller(c(CTR = "Vehicle", LEN = "Lenvatinib"))) +
  scale_color_manual(values = cluster_palette, drop = FALSE) +
  native_coord() +
  labs(
    title = "B  De novo TAN clusters",
    subtitle = "Same independently calculated query-native Neu UMAP",
    x = "UMAP 1", y = "UMAP 2", color = "Cluster"
  ) + umap_theme

p_ng <- ggplot(per_cell, aes(UMAP_1, UMAP_2, color = Ng_stage)) +
  geom_point(size = 0.20, alpha = 0.80, stroke = 0) +
  scale_color_manual(values = ng_palette, drop = FALSE) +
  native_coord() +
  labs(title = "C  Ng maturation labels", x = "UMAP 1", y = "UMAP 2", color = "Ng stage") +
  umap_theme
p_xue <- ggplot(per_cell, aes(UMAP_1, UMAP_2, color = Xue_state)) +
  geom_point(size = 0.20, alpha = 0.80, stroke = 0) +
  scale_color_manual(values = xue_palette, drop = FALSE) +
  native_coord() +
  labs(
    title = "Transferred Xue tumour-mNeu states",
    x = "UMAP 1", y = "UMAP 2", color = "Xue state"
  ) + umap_theme
p_reference_labels <- p_ng + p_xue +
  plot_layout(ncol = 2) +
  plot_annotation(
    subtitle = paste0(
      "Discrete reference labels on identical query-native coordinates; ",
      "Xue labels use pTMC/pTMK tumour states 07-12"
    )
  )

ridge_data <- per_cell |>
  mutate(
    Neu_cluster_ridge = factor(
      as.character(Neu_cluster), levels = rev(cluster_display_levels)
    )
  )
p_maturation <- ggplot(
  ridge_data,
  aes(Xue_core_maturation, Neu_cluster_ridge, fill = Neu_cluster_ridge)
) +
  geom_hline(
    yintercept = seq_along(cluster_display_levels),
    color = "black", linewidth = 0.35
  ) +
  ggridges::geom_density_ridges(
    scale = 1.18, alpha = 0.95, color = "black", linewidth = 0.35,
    rel_min_height = 0.01
  ) +
  scale_fill_manual(values = cluster_palette, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, 1), breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "D  Neutrophil maturation signature",
    subtitle = "Xue core-state marker/UCell axis; C1 at top and C5 at bottom",
    x = "Neutrophil maturation signature score", y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none", plot.title = element_text(face = "bold")
  )

p_composition <- cluster_composition |>
  filter(group %in% c("CTR", "LEN")) |>
  mutate(group = factor(as.character(group), levels = c("CTR", "LEN"))) |>
  ggplot(aes(group, percent_within_captured_Neu, fill = Neu_cluster)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = cluster_palette, drop = FALSE) +
  scale_x_discrete(labels = c(CTR = "Veh", LEN = "Len")) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "E  Recovered TAN cluster composition",
    subtitle = "Descriptive; not intratumoral abundance",
    x = NULL, y = "% within recovered neutrophils", fill = "Cluster"
  ) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))
p_score_composition <- p_maturation + p_composition +
  plot_layout(widths = c(1.65, 1))

p_neu_combined <- p_marker_dotplot / p_cluster_ctr_len /
  p_reference_labels / p_score_composition +
  plot_layout(heights = c(1.20, 1, 1, 1)) +
  plot_annotation(
    title = "Refined characterization of lenvatinib-associated TAN states",
    caption = paste0(
      "One pooled 10x library per condition. Compositions and score distributions ",
      "are descriptive and do not constitute treatment-level statistical tests."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  )

figure_files <- c(
  save_plot_both(p_marker_dotplot, "Fig8A_de_novo_marker_DotPlot", 14.5, 5.8),
  save_plot_both(p_cluster_ctr_len, "Fig8B_CTR_LEN_de_novo_query_native_UMAP", 11.5, 5.5),
  save_plot_both(p_reference_labels, "Fig8C_Ng_Xue_labels_query_native_UMAP", 13.0, 5.7),
  save_plot_both(p_score_composition, "Fig8D_E_Xue_core_maturation_composition", 12.5, 5.4),
  save_plot_both(p_neu_combined, "Fig8A_E_refined_TAN_characterization", 17.0, 20.0)
)

# %% 10 - save complete audit outputs and final object

stage("Writing tables, mapped object, parameters, and session information.")
table_files <- c(
  input_manifest_file,
  write_csv(balanced_hvg$ranking, "balanced_query_HVG_ranking"),
  write_csv(balanced_hvg$within_group, "balanced_query_HVGs_by_library"),
  write_csv(xue_reference_audit, "Xue_all_mNeu_sample_tissue_state_audit"),
  write_csv(xue_reference_counts, "Xue_tumour_reference_counts"),
  write_csv(xue_numeric_crosswalk, "Xue_numeric_state_crosswalk"),
  write_csv(mapping_qc_overall, "Xue_mapping_QC_overall"),
  write_csv(mapping_qc_by_group, "Xue_mapping_QC_by_group"),
  write_csv(xue_composition, "Xue_state_composition_by_group"),
  write_csv(ng_xue_crosswalk, "Ng_Xue_state_crosswalk"),
  write_csv(xue_core_marker_manifest, "Xue_core_marker_manifest_all"),
  write_csv(xue_core_coverage, "Xue_core_marker_coverage"),
  write_csv(cluster_relabel, "de_novo_cluster_relabel_crosswalk"),
  write_csv(cluster_sizes, "de_novo_cluster_sizes"),
  write_csv(cluster_composition, "de_novo_cluster_composition_by_group"),
  write_csv(maturation_summary, "Xue_core_maturation_summary_by_cluster"),
  write_csv(marker_manifest, "de_novo_selected_marker_manifest"),
  write_csv(marker_dotplot_data, "de_novo_marker_DotPlot_data"),
  write_csv(cluster_markers, "de_novo_FindAllMarkers_all"),
  write_csv(cluster_markers_top20, "de_novo_FindAllMarkers_top20_by_cluster"),
  write_csv(per_cell, "per_cell_metadata")
)

parameter_manifest <- tibble(
  section = c(
    rep("query_native_UMAP", 7), rep("Xue_mapping", 11),
    rep("Xue_core_maturation", 7), rep("de_novo_clustering", 8),
    rep("interpretation", 3)
  ),
  parameter = c(
    "HVG_strategy", "HVG_per_library", "HVG_target", "PCA_npcs",
    "PCA_dims", "UMAP", "seed",
    "reference_source", "reference_celltype", "reference_samples",
    "reference_tissue", "reference_states", "reference_cells", "HVG_target",
    "anchor_reduction", "PCA_dims", "query_display_geometry", "seed",
    "marker_source_cells", "early_states", "mature_like_states",
    "marker_filter", "markers_per_state", "scoring", "display_transform",
    "assay", "PCA_npcs", "SNN_dims", "k_param", "algorithm", "resolution",
    "seed", "relabel",
    "pooled_library_design", "composition_inference", "cell_level_treatment_tests"
  ),
  value = c(
    "within-library rank consensus; no treatment/batch integration", "4000",
    "3000", "50", "1:30", "n.neighbors=50;min.dist=0.35;metric=cosine", "123",
    xue_file, "mNeu", "pTMC and pTMK", "Tumor", "mNeu states 07-12", "8297",
    "2000 reference HVGs shared with query", "pcaproject", "1:30",
    "independently calculated query-native Neu UMAP", as.character(mapping_seed),
    "all 17,780 Xue mouse mNeu cells", "01_Ngp and 02_Mmp8",
    "03_Pabpc1,05_Gm2a,06_Marco,07_Actg1",
    "presto one-vs-rest: padj<0.05;logFC>=0.25;pct_in>=0.10;AUC>=0.55",
    "top 30 per state", "UCell", "rank within score then mature-minus-early; rescale 0-1",
    "integrated", "30", "1:20", "50", "Louvain algorithm 1", "0.3",
    as.character(cluster_seed), "frozen core-state-informed C5,C3,C4,C1,C2 -> C1-C5",
    "one pooled 10x library per condition", "descriptive only", "none"
  )
)
parameter_file <- write_csv(parameter_manifest, "parameter_manifest")
table_files <- c(table_files, parameter_file)

object_file <- save_rds(
  query,
  "Len_Neu_Xue_tumour_mapping_Xue_core_score_de_novo_res0p3",
  compress = TRUE
)
session_file <- file.path(dirs$logs, paste0(prefix, "__sessionInfo.txt"))
capture.output(sessionInfo(), file = session_file)
command_file <- file.path(dirs$logs, paste0(prefix, "__command.txt"))
writeLines(c(
  paste(commandArgs(), collapse = " "),
  paste0("run_dir=", normalizePath(run_dir, winslash = "/", mustWork = FALSE))
), command_file)

core_outputs <- unique(c(
  table_files, figure_files, object_file, session_file, command_file
))
if (any(!file.exists(core_outputs)) || any(file.info(core_outputs)$size <= 0)) {
  stop("At least one expected core output is missing or empty.")
}

validation <- tibble(
  gate = c(
    "query_has_14111_Neu_cells", "six_pooled_libraries_present",
    "query_native_UMAP_complete", "Xue_all_mNeu_17780",
    "Xue_tumour_reference_8297", "Xue_states_07_12_present",
    "Xue_label_score_agreement", "Xue_core_score_complete",
    "five_de_novo_clusters", "C1_top_C5_bottom_display_contract",
    "no_cell_level_treatment_tests", "core_outputs_nonempty"
  ),
  passed = c(
    ncol(query) == 14111L,
    setequal(unique(as.character(query$group)), group_levels),
    nrow(query_native_umap_after_clustering) == 14111L &&
      all(is.finite(query_native_umap_after_clustering)),
    xue_all_mneu_n == 17780L,
    ncol(xue_ref) == 8297L,
    setequal(unique(as.character(query$Xue_state)), xue_tumour_state_levels),
    label_score_agreement >= 0.999,
    all(is.finite(query$Xue_core_maturation)),
    nrow(cluster_sizes) == 5L && all(cluster_sizes$n_cells > 0L),
    identical(levels(ridge_data$Neu_cluster_ridge), rev(cluster_display_levels)),
    TRUE,
    all(file.exists(core_outputs)) && all(file.info(core_outputs)$size > 0)
  ),
  severity = "ERROR"
)
validation_file <- write_csv(validation, "validation")
if (any(!validation$passed)) {
  stop("Validation failed: ", paste(validation$gate[!validation$passed], collapse = ", "))
}

writeLines(
  c(
    "status=COMPLETED",
    paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("query_cells=", ncol(query)),
    paste0("Xue_all_mNeu_cells=", xue_all_mneu_n),
    paste0("Xue_tumour_reference_cells=", ncol(xue_ref)),
    paste0("de_novo_clusters=", nrow(cluster_sizes)),
    "treatment_level_inference=none_one_pooled_library_per_condition",
    paste0("validation_file=", validation_file)
  ),
  status_file
)

stage("Validation passed; finalizing output manifest: ", run_dir)
print(validation)
close_log()

manifest_inputs <- unique(c(
  core_outputs, validation_file, status_file, log_file
))
if (any(!file.exists(manifest_inputs)) || any(file.info(manifest_inputs)$size <= 0)) {
  stop("At least one final output is missing or empty.")
}
output_manifest <- tibble(
  path = normalizePath(manifest_inputs, winslash = "/", mustWork = TRUE),
  size_bytes = file.info(manifest_inputs)$size,
  md5 = unname(tools::md5sum(manifest_inputs))
)
output_manifest_file <- write_csv(output_manifest, "output_manifest")
if (!file.exists(output_manifest_file) || file.info(output_manifest_file)$size <= 0) {
  stop("Output manifest is missing or empty.")
}

run_complete <- TRUE
cat("STEP6_PASS\n")
cat("output_dir\t", normalizePath(run_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")

}

main()

