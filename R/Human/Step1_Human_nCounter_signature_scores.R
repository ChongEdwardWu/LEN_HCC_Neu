#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste(
    "Usage: Rscript Step1_Human_nCounter_signature_scores.R",
    "<normalized_nCounter.xlsx> <sample_metadata.csv> <output_dir> [sheet]"
  )
  if (length(args) == 1L && args %in% c("--help", "-h")) {
    cat(usage, "\nMetadata columns: sample_id, patient_id, timepoint (baseline/week4/PD).\n")
    return(invisible(NULL))
  }
  if (!length(args) %in% 3:4) stop(usage, call. = FALSE)
  output_dir <- args[[3]]
  sheet <- if (length(args) == 4L) args[[4]] else "NC_9785_N1_Normalization"
  for (pkg in c("readxl", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package: ", pkg)
  }
  if (file.exists(output_dir) && (!dir.exists(output_dir) ||
      length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))) {
    stop("Output directory must be new or empty; no files were changed.")
  }
  gene_sets <- list(
    B_cell = c("CD19", "MS4A1"), CD3_T = c("CD3D", "CD3E", "CD3G"),
    CD4_T = "CD4", CD8_T = c("CD8A", "CD8B"), Macrophage = c("CD68", "CSF1R"),
    Monocyte = c("CD14", "CCR2"), NK = c("KLRD1", "KLRK1", "NCR1", "NKG7", "GNLY"),
    Neutrophil = c("CXCR2", "CSF3R"), cDC = "ITGAX", pDC = "IL3RA"
  )

  raw <- suppressMessages(readxl::read_xlsx(
    args[[1]], sheet = sheet, col_names = FALSE, col_types = "text"
  ))
  if (nrow(raw) < 5L || ncol(raw) < 8L) stop("Unexpected nCounter workbook layout.")
  header <- trimws(as.character(unlist(raw[1, ], use.names = FALSE)))
  sample_cols <- seq.int(8L, ncol(raw))
  ids <- trimws(as.character(unlist(raw[3, sample_cols], use.names = FALSE)))
  sample_cols <- sample_cols[!is.na(ids) & grepl("^FL", ids)]
  sample_ids <- trimws(as.character(unlist(raw[3, sample_cols], use.names = FALSE)))
  if (!length(sample_ids) || anyDuplicated(sample_ids)) stop("Missing or duplicate sample IDs.")
  # Both supplied workflows must refer to the same set of sample columns.
  rcc_cols <- which(!is.na(header) & grepl("\\.RCC$", header))
  if (!setequal(rcc_cols, sample_cols)) {
    stop("RCC-header and FL-ID sample selections differ. Verify the workbook sample columns.")
  }
  if (is.na(header[[1]]) || header[[1]] != "Probe Name") stop("Column 1 must be Probe Name.")
  genes <- trimws(as.character(raw[[1]][5:nrow(raw)]))
  markers <- unique(unlist(gene_sets, use.names = FALSE))
  if (length(setdiff(markers, genes))) {
    stop("Missing marker genes: ", paste(setdiff(markers, genes), collapse = ", "))
  }
  if (anyDuplicated(genes[genes %in% markers])) stop("Duplicate marker-gene rows.")
  marker_rows <- match(markers, genes) + 4L
  expr <- as.matrix(raw[marker_rows, sample_cols, drop = FALSE])
  suppressWarnings(storage.mode(expr) <- "numeric")
  if (any(!is.finite(expr)) || any(expr < 0)) stop("Marker expression must be finite and nonnegative.")
  dimnames(expr) <- list(markers, sample_ids)

  meta <- utils::read.csv(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sample_id", "patient_id", "timepoint")
  if (anyDuplicated(names(meta)) || !all(required %in% names(meta))) {
    stop("Metadata requires unique columns: sample_id, patient_id, timepoint.")
  }
  meta <- meta[, required, drop = FALSE]
  meta[] <- lapply(meta, function(x) trimws(as.character(x)))
  if (anyNA(meta) || any(as.matrix(meta) == "") || anyDuplicated(meta$sample_id)) {
    stop("Metadata IDs and timepoints must be nonempty, with unique sample_id values.")
  }
  if (!setequal(meta$sample_id, sample_ids)) stop("Metadata and nCounter sample IDs must match exactly.")
  if (!all(meta$timepoint %in% c("baseline", "week4", "PD"))) stop("Unknown metadata timepoint.")
  if (anyDuplicated(meta[c("patient_id", "timepoint")])) stop("Duplicate patient/timepoint records.")
  meta <- meta[match(sample_ids, meta$sample_id), , drop = FALSE]

  expr_log <- log2(expr + 1)
  scores <- vapply(gene_sets, function(g) colMeans(expr_log[g, , drop = FALSE]),
                   numeric(ncol(expr)))
  baseline <- meta$timepoint == "baseline"
  if (sum(baseline) < 2L) stop("At least two baseline samples are required.")
  mu <- colMeans(scores[baseline, , drop = FALSE])
  sigma <- apply(scores[baseline, , drop = FALSE], 2, stats::sd)
  if (any(!is.finite(sigma) | sigma <= 0)) stop("Every baseline signature SD must be positive.")
  standardized <- sweep(sweep(scores, 2, mu, "-"), 2, sigma, "/")
  colnames(standardized) <- paste0(names(gene_sets), "_z")
  score_df <- cbind(meta, as.data.frame(scores), as.data.frame(standardized))

  threshold_col <- which(header == "% Samples above Threshold")
  if (length(threshold_col) > 1L) stop("Duplicate threshold-percentage columns.")
  if (length(threshold_col)) {
    pct <- suppressWarnings(as.numeric(raw[[threshold_col]][marker_rows]))
    if (any(!is.finite(pct)) || any(pct < 0 | pct > 100)) {
      stop("Marker threshold percentages must be present and within [0, 100].")
    }
    fraction <- pct / 100
    size_label <- "Fraction above threshold"
    detection_source <- "Workbook percentage above threshold"
  } else {
    fraction <- rowMeans(expr > 0)
    size_label <- "Fraction with expression > 0"
    detection_source <- "Fraction of samples with expression > 0"
  }
  marker_summary <- data.frame(
    gene = markers, mean_expr = rowMeans(expr),
    mean_log2 = log2(rowMeans(expr) + 1), pct_use = fraction,
    detection_source = detection_source, stringsAsFactors = FALSE
  )
  plot_df <- data.frame(
    signature = rep(names(gene_sets), lengths(gene_sets)),
    gene = unlist(gene_sets, use.names = FALSE), stringsAsFactors = FALSE
  )
  plot_df <- cbind(plot_df, marker_summary[match(plot_df$gene, marker_summary$gene), -1])
  plot_df$signature <- factor(plot_df$signature, levels = names(gene_sets))
  plot_df$gene <- factor(plot_df$gene, levels = marker_summary$gene[order(-marker_summary$mean_log2)])
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = gene, y = signature)) +
    ggplot2::geom_point(ggplot2::aes(size = pct_use, color = mean_log2), alpha = 0.95) +
    ggplot2::scale_size_continuous(range = c(1.5, 7)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "Immune signature definitions with expression overlay (all samples)",
      x = "Marker gene", y = "Immune signature",
      color = "log2(mean expression + 1)", size = size_label
    )
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) stop("Cannot create output directory.")
  utils::write.csv(score_df, file.path(output_dir, "signature_scores_and_z_by_sample.csv"), row.names = FALSE)
  utils::write.csv(data.frame(signature = names(mu), baseline_mean = mu, baseline_sd = sigma),
                   file.path(output_dir, "baseline_signature_reference.csv"), row.names = FALSE)
  utils::write.csv(plot_df, file.path(output_dir, "Signature_marker_expression_summary_ALLsamples.csv"), row.names = FALSE)
  ggplot2::ggsave(file.path(output_dir, "Signature_definition_dot_with_expression_ALLsamples.pdf"),
                 p, width = 10.5, height = 5.2)
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  message("Scored ", nrow(meta), " samples using ", sum(baseline), " baseline samples.")
  invisible(score_df)
}

if (sys.nframe() == 0L) main()
