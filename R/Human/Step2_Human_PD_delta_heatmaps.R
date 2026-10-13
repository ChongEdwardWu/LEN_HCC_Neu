#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste(
    "Usage: Rscript Step2_Human_PD_delta_heatmaps.R",
    "<patient_annotation.xlsx> <signature_scores_and_z_by_sample.csv> <output_dir> [sheet]"
  )
  if (length(args) == 1L && args %in% c("--help", "-h")) {
    cat(usage, "\nClinical columns: FL_ID, mRECIST. Default sheet: Patient annotation.\n")
    return(invisible(NULL))
  }
  if (!length(args) %in% 3:4) stop(usage, call. = FALSE)
  output_dir <- args[[3]]
  sheet <- if (length(args) == 4L) args[[4]] else "Patient annotation"
  for (pkg in c("readxl", "pheatmap")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package: ", pkg)
  }
  if (file.exists(output_dir) && (!dir.exists(output_dir) ||
      length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))) {
    stop("Output directory must be new or empty; no files were changed.")
  }
  clin <- as.data.frame(readxl::read_xlsx(args[[1]], sheet = sheet))
  if (!all(c("FL_ID", "mRECIST") %in% names(clin))) stop("Clinical columns FL_ID and mRECIST are required.")
  patient_id <- trimws(as.character(clin$FL_ID))
  if (anyNA(patient_id) || any(!nzchar(patient_id)) || anyDuplicated(patient_id)) stop("Clinical patient IDs must be unique and nonempty.")
  response <- toupper(trimws(as.character(clin$mRECIST)))
  group <- ifelse(response == "PD", "PD", ifelse(response %in% c("CR", "PR", "SD"), "non-PD", NA_character_))
  names(group) <- patient_id

  scores <- utils::read.csv(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
  for (alias in c("PMN_z", "PMN_MDSC_z")) {
    if (alias %in% names(scores)) {
      if ("Neutrophil_z" %in% names(scores)) stop("Conflicting neutrophil signature aliases.")
      names(scores)[names(scores) == alias] <- "Neutrophil_z"
    }
  }
  signatures <- c("B_cell", "CD3_T", "CD4_T", "CD8_T", "Macrophage", "Monocyte", "NK", "Neutrophil", "cDC", "pDC")
  z_cols <- paste0(signatures, "_z")
  if (anyDuplicated(names(scores)) || !all(c("sample_id", "patient_id", "timepoint", z_cols) %in% names(scores))) {
    stop("Scores must contain sample_id, patient_id, timepoint, and all ten signature Z columns.")
  }
  if (anyNA(scores[c("sample_id", "patient_id", "timepoint")]) ||
      any(as.matrix(scores[c("sample_id", "patient_id", "timepoint")]) == "") ||
      anyDuplicated(scores$sample_id) || anyDuplicated(scores[c("patient_id", "timepoint")])) {
    stop("Missing or duplicate sample/patient/timepoint identifiers.")
  }
  if (!all(scores$timepoint %in% c("baseline", "week4", "PD"))) stop("Unknown timepoint in score table.")
  for (column in z_cols) {
    scores[[column]] <- suppressWarnings(as.numeric(scores[[column]]))
    if (any(!is.finite(scores[[column]]))) stop("Nonfinite Z scores in ", column)
  }
  # PD denotes the biopsy timepoint here; mRECIST defines the separate response group.
  baseline <- scores[scores$timepoint == "baseline", , drop = FALSE]
  followup <- scores[scores$timepoint == "PD", , drop = FALSE]
  paired_ids <- intersect(baseline$patient_id, followup$patient_id)
  paired_ids <- paired_ids[!is.na(group[paired_ids])]
  if (!length(paired_ids)) stop("No complete baseline/PD pairs with an mRECIST group.")
  baseline <- baseline[match(paired_ids, baseline$patient_id), , drop = FALSE]
  followup <- followup[match(paired_ids, followup$patient_id), , drop = FALSE]
  paired_group <- unname(group[paired_ids])
  if (!all(c("non-PD", "PD") %in% paired_group)) stop("Both mRECIST response groups are required.")
  delta_pd <- do.call(rbind, lapply(seq_along(signatures), function(i) data.frame(
    patient_id = paired_ids, group = paired_group, signature = signatures[[i]],
    baseline = baseline[[z_cols[[i]]]], PD = followup[[z_cols[[i]]]],
    delta = followup[[z_cols[[i]]]] - baseline[[z_cols[[i]]]],
    group2 = ifelse(paired_group == "non-PD", "nonPD", "PD"), stringsAsFactors = FALSE
  )))
  group_mean <- do.call(rbind, lapply(signatures, function(s) {
    do.call(rbind, lapply(c("PD", "nonPD"), function(g) {
      x <- delta_pd$delta[delta_pd$signature == s & delta_pd$group2 == g]
      data.frame(signature = s, group2 = g, n = length(x),
                 mean_delta = mean(x), median_delta = stats::median(x))
    }))
  }))
  wide <- data.frame(signature = signatures)
  for (metric in c("n", "mean_delta", "median_delta")) {
    for (g in c("PD", "nonPD")) {
      rows <- group_mean[group_mean$group2 == g, , drop = FALSE]
      wide[[paste0(metric, "__", g)]] <- rows[[metric]][match(wide$signature, rows$signature)]
    }
  }
  patient_order <- order(paired_group, paired_ids)
  patient_mat <- t(as.matrix(followup[z_cols]) - as.matrix(baseline[z_cols]))
  dimnames(patient_mat) <- list(signatures, paired_ids)
  patient_mat <- patient_mat[signatures, patient_order, drop = FALSE]
  annotation <- data.frame(group = paired_group[patient_order], row.names = paired_ids[patient_order])
  mean_mat <- as.matrix(wide[c("mean_delta__PD", "mean_delta__nonPD")])
  dimnames(mean_mat) <- list(wide$signature, c("PD", "non-PD"))
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) stop("Cannot create output directory.")
  utils::write.csv(delta_pd, file.path(output_dir, "PD_delta_long_by_patient_signature.csv"), row.names = FALSE)
  utils::write.csv(group_mean, file.path(output_dir, "PD_group_mean_median_delta_long.csv"), row.names = FALSE)
  utils::write.csv(wide, file.path(output_dir, "PD_group_mean_median_delta_wide.csv"), row.names = FALSE)
  grDevices::png(file.path(output_dir, "PD_delta_heatmap.png"), width = 1200, height = 800, res = 150)
  tryCatch(pheatmap::pheatmap(patient_mat, cluster_rows = TRUE, cluster_cols = FALSE,
                            annotation_col = annotation, main = "PD biopsy - baseline: signature Delta Z"),
           finally = grDevices::dev.off())
  grDevices::png(file.path(output_dir, "PD_group_mean_heatmap.png"), width = 900, height = 700, res = 150)
  tryCatch(pheatmap::pheatmap(mean_mat, cluster_rows = TRUE, cluster_cols = FALSE,
                            main = "Mean Delta Z by mRECIST response group"),
           finally = grDevices::dev.off())
  pheatmap::pheatmap(mean_mat, cluster_rows = TRUE, cluster_cols = FALSE,
                    main = "Mean Delta Z by mRECIST response group",
                    filename = file.path(output_dir, "PD_group_mean_heatmap.pdf"), width = 6, height = 6)
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  message("Paired patients by mRECIST response group:")
  print(table(paired_group))
  invisible(list(delta = delta_pd, summary = group_mean))
}

if (sys.nframe() == 0L) main()
