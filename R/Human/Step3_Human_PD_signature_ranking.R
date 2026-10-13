usage <- function() {
  cat(paste0(
    "Usage: Rscript Step3_Human_PD_signature_ranking.R ",
    "<STable_PD_delta_stats_ALLsignatures.csv> <output_dir>\n\n",
    "Descriptively ranks supplied PD versus non-PD signature statistics.\n",
    "The output directory must not exist or must be empty.\n",
    "Required package: ggplot2. Confidence intervals are not recalculated.\n"
  ))
}

z <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

read_signature_stats <- function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    stop("Input CSV does not exist or is not a file: ", path, call. = FALSE)
  }
  stats <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "signature", "diff_median", "prop_up_PD", "prop_up_nonPD",
    "r_rank_biserial", "ci_med_low", "ci_med_high"
  )
  if (anyDuplicated(names(stats))) {
    stop("Input CSV has duplicate column names.", call. = FALSE)
  }
  missing <- setdiff(required, names(stats))
  if (length(missing)) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!nrow(stats)) stop("Input CSV has no data rows.", call. = FALSE)

  stats$signature <- trimws(as.character(stats$signature))
  if (anyNA(stats$signature) || any(!nzchar(stats$signature))) {
    stop("Signature names must be nonmissing and nonempty.", call. = FALSE)
  }
  if (anyDuplicated(stats$signature)) {
    stop("Each signature must occur exactly once.", call. = FALSE)
  }
  aliases <- stats$signature %in% c("PMN", "PMN_MDSC")
  if (any(aliases)) {
    message("Renaming legacy signature aliases to Neutrophil: ",
            paste(stats$signature[aliases], collapse = ", "))
    stats$signature[aliases] <- "Neutrophil"
  }
  if (anyDuplicated(stats$signature)) {
    stop("Signature alias collision after mapping PMN/PMN_MDSC to Neutrophil.",
         call. = FALSE)
  }

  for (column in setdiff(required, "signature")) {
    value <- suppressWarnings(as.numeric(stats[[column]]))
    if (any(!is.finite(value))) {
      stop("Column must contain finite numeric values: ", column, call. = FALSE)
    }
    stats[[column]] <- value
  }
  for (column in c("prop_up_PD", "prop_up_nonPD")) {
    if (any(stats[[column]] < 0 | stats[[column]] > 1)) {
      stop("Proportions must be within [0, 1]: ", column, call. = FALSE)
    }
  }
  if (any(abs(stats$r_rank_biserial) > 1)) {
    stop("Absolute rank-biserial correlations must not exceed 1.", call. = FALSE)
  }
  if (any(stats$ci_med_low > stats$ci_med_high)) {
    stop("Every ci_med_low must be less than or equal to ci_med_high.", call. = FALSE)
  }
  stats
}

rank_signatures <- function(stats) {
  # Standardization is across the supplied signatures, not across patients.
  stats$diff_prop_up <- stats$prop_up_PD - stats$prop_up_nonPD
  stats$r_use <- sign(stats$diff_median) * abs(stats$r_rank_biserial)
  stats$z_diff_median <- z(stats$diff_median)
  stats$z_diff_prop_up <- z(stats$diff_prop_up)
  stats$z_r_use <- z(stats$r_use)
  stats$score <- 0.5 * stats$z_diff_median +
    0.2 * stats$z_diff_prop_up + 0.3 * stats$z_r_use
  if (any(!is.finite(stats$score))) {
    stop("Composite scores are nonfinite; check the scale of the input values.",
         call. = FALSE)
  }
  stats <- stats[order(-stats$score), , drop = FALSE]
  rownames(stats) <- NULL
  stats
}

signature_lollipop <- function(stats, myeloid = FALSE) {
  stats$signature <- stats::reorder(factor(stats$signature), stats$score)
  title <- if (myeloid) {
    "Myeloid-focused ranking (PD vs non-PD)"
  } else {
    "Ranking of immune signatures by composite score (PD vs non-PD)"
  }
  subtitle <- if (myeloid) {
    "Prioritization by composite score\nThick line = supplied bootstrap 95% CI"
  } else {
    "X = median difference of Delta Z (PD - non-PD)\nThick line = supplied bootstrap 95% CI"
  }
  ggplot2::ggplot(stats, ggplot2::aes(x = diff_median, y = signature)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = diff_median, yend = signature), linewidth = 0.6
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = ci_med_low, xend = ci_med_high, yend = signature),
      linewidth = 1.0
    ) +
    ggplot2::geom_point(size = if (myeloid) 2.6 else 2.4) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title, subtitle = subtitle,
      x = "Median difference of Delta Z (PD - non-PD)", y = NULL
    )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 1L && args %in% c("--help", "-h")) {
    usage()
    return(invisible(NULL))
  }
  if (length(args) != 2L) {
    usage()
    stop("Exactly two positional arguments are required.", call. = FALSE)
  }
  input_file <- args[[1L]]
  output_dir <- args[[2L]]
  if (!nzchar(output_dir)) stop("Output directory must not be empty.", call. = FALSE)
  if (file.exists(output_dir) && !dir.exists(output_dir)) {
    stop("Output path exists and is not a directory.", call. = FALSE)
  }
  if (dir.exists(output_dir) &&
      length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory is not empty; no files were changed.", call. = FALSE)
  }

  rank_tbl <- rank_signatures(read_signature_stats(input_file))
  myeloid <- c("Neutrophil", "Macrophage", "Monocyte", "cDC", "pDC")
  plot_my <- rank_tbl[rank_tbl$signature %in% myeloid, , drop = FALSE]
  if (!nrow(plot_my)) {
    stop("Input contains no recognized myeloid signatures.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required. Install it before running this script.",
         call. = FALSE)
  }
  p_all <- signature_lollipop(rank_tbl)
  p_my <- signature_lollipop(plot_my, myeloid = TRUE)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop("Could not create output directory: ", output_dir, call. = FALSE)
  }

  utils::write.csv(rank_tbl, file.path(output_dir, "PD_signature_rank_composite_score.csv"),
                   row.names = FALSE)
  ggplot2::ggsave(file.path(output_dir, "SFig_rank_lollipop_ALL.pdf"),
                 p_all, width = 8, height = 5.5)
  ggplot2::ggsave(file.path(output_dir, "SFig_rank_lollipop_MYeloid.pdf"),
                 p_my, width = 6.5, height = 4.2)
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  message("Descriptive ranking completed for ", nrow(rank_tbl), " signatures. Outputs: ",
          normalizePath(output_dir, winslash = "/", mustWork = TRUE))
  invisible(rank_tbl)
}

if (sys.nframe() == 0L) main()
