# Neutrophil-to-CD8Tc nearest distances, calculated within each QupathID.
# Usage: Rscript Distance_and_phenotype_analysis.R cells.rds output_dir

normalize_distance <- function(x) {
  result <- rep(NA_real_, length(x))
  valid <- is.finite(x)
  if (any(valid) && max(x[valid]) > 0) {
    result[valid] <- 100 * x[valid] / max(x[valid])
  }
  result
}

calculate_neutrophil_distances <- function(csd,
    distance_fun = phenoptr::find_nearest_distance) {
  if (!is.data.frame(csd) || nrow(csd) == 0L) {
    stop("The RDS input must contain a non-empty cell data frame.")
  }
  csd <- as.data.frame(csd, check.names = FALSE)
  if (!"Cell ID" %in% names(csd) && "CellID" %in% names(csd)) {
    names(csd)[names(csd) == "CellID"] <- "Cell ID"
  }
  required <- c("QupathID", "SampleName", "Cell ID", "Cell X Position",
                "Cell Y Position", "Phenotype")
  missing <- setdiff(required, names(csd))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  if (anyNA(csd[, required])) stop("Required cell-data columns must be complete.")
  for (column in c("Cell X Position", "Cell Y Position")) {
    if (!is.numeric(csd[[column]]) || any(!is.finite(csd[[column]]))) {
      stop(column, " must contain finite numeric coordinates.")
    }
  }
  if (!all(c("Neutrophil", "CD8Tc") %in% csd$Phenotype)) {
    stop("The input must include the phenotype labels Neutrophil and CD8Tc.")
  }

  # Keep each sample separate for both distance calculation and normalization.
  sample_results <- lapply(unique(csd$QupathID), function(sample_id) {
    cells <- csd[csd$QupathID == sample_id, , drop = FALSE]
    if (anyDuplicated(cells[["Cell ID"]])) {
      stop("Cell IDs must be unique within QupathID: ", sample_id)
    }
    # phenoptr uses this column to identify a single spatial field.
    cells[["Sample Name"]] <- as.character(sample_id)
    distances <- distance_fun(cells, phenotypes = list(CD8Tc = "CD8Tc"))
    distance_column <- intersect(c("Distance to CD8Tc", "Distance.to.CD8Tc"),
                                 names(distances))
    if (length(distance_column) != 1L || nrow(distances) != nrow(cells)) {
      stop("Unexpected nearest-distance output for QupathID: ", sample_id)
    }
    keep <- cells$Phenotype == "Neutrophil"
    result <- cells[keep, required, drop = FALSE]
    result$Distance.to.CD8Tc <- distances[[distance_column]][keep]
    result$normalized_CD8Tc <- normalize_distance(result$Distance.to.CD8Tc)
    result
  })
  result <- do.call(rbind, sample_results)
  rownames(result) <- NULL
  result
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L || any(args %in% c("-h", "--help"))) {
    stop("Usage: Rscript Distance_and_phenotype_analysis.R cells.rds output_dir",
         call. = FALSE)
  }
  packages <- c("phenoptr", "rtree")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
  if (!file.exists(args[1])) stop("Input RDS file does not exist: ", args[1])
  output_dir <- args[2]
  if (dir.exists(output_dir) &&
      length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Use a new or empty output directory: ", output_dir)
  }
  result <- calculate_neutrophil_distances(readRDS(args[1]))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(result, file.path(output_dir, "neutrophil_CD8Tc_distances.csv"),
            row.names = FALSE, na = "NA")
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  message("Saved distances for ", nrow(result), " neutrophils to ", output_dir)
}

if (sys.nframe() == 0L) main()
