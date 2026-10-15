# TCGA-LIHC overall survival by CIBERSORT neutrophils and PERK-UPR ssGSEA.
# Usage: Rscript KM_survival.R expression_log2_tpm.tsv survival.tsv output_dir
# Expression: Ensembl IDs in column 1; TCGA sample IDs in subsequent columns;
# values on the log2(TPM + 1) scale used by the supplied analysis.
# Survival: sample, OS (0 = censored; 1 = event), OS.time (days).

prepare_survival_data <- function(score_df, clinical) {
  required <- c("sample", "OS", "OS.time")
  missing <- setdiff(required, names(clinical))
  if (length(missing)) stop("Missing survival columns: ", paste(missing, collapse = ", "))
  if (anyNA(clinical$sample) || anyDuplicated(clinical$sample)) {
    stop("Survival sample IDs must be complete and unique.")
  }
  result <- merge(score_df, clinical[, required, drop = FALSE], by = "sample", all = FALSE)
  if (!nrow(result)) stop("Expression and survival sample IDs must match.")
  for (column in c("OS", "OS.time", "Neu_sig", "PERK_UPR")) {
    if (!is.numeric(result[[column]]) || any(!is.finite(result[[column]]))) {
      stop(column, " must contain complete finite numeric values in matched samples.")
    }
  }
  if (any(!result$OS %in% c(0, 1)) || any(result$OS.time < 0)) {
    stop("OS must be 0 or 1, and OS.time must be non-negative days.")
  }
  result$OS.day <- result$OS.time
  result$OS.time <- result$OS.time / 365 * 12
  result$OS <- ifelse(result$OS.time > 72, 0, result$OS)
  result$OS.time <- pmin(result$OS.time, 72)
  result
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L || any(args %in% c("-h", "--help"))) {
    stop("Usage: Rscript KM_survival.R expression_log2_tpm.tsv survival.tsv output_dir",
         call. = FALSE)
  }
  packages <- c("data.table", "AnnotationDbi", "org.Hs.eg.db", "IOBR",
                "msigdbr", "GSVA", "survival", "survminer", "maxstat", "ggplot2")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages(library(survival))
  if (!all(file.exists(args[1:2]))) stop("Both input TSV files must exist.")
  output_dir <- args[3]
  if (dir.exists(output_dir) &&
      length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Use a new or empty output directory: ", output_dir)
  }
  set.seed(20261015)

  # 1. Read log2(TPM + 1), clean Ensembl IDs, and retain the highest-mean row.
  tpm_dt <- data.table::fread(args[1], sep = "\t", header = TRUE, data.table = FALSE)
  if (ncol(tpm_dt) < 2L || nrow(tpm_dt) == 0L) stop("Expression table is empty.")
  tpm <- as.matrix(tpm_dt[, -1, drop = FALSE])
  rownames(tpm) <- as.character(tpm_dt[[1]])
  if (!is.numeric(tpm) || any(!is.finite(tpm)) || any(tpm < 0)) {
    stop("Expression values must be finite, non-negative log2(TPM + 1) values.")
  }
  if (anyNA(rownames(tpm)) || any(!nzchar(rownames(tpm))) ||
      anyDuplicated(colnames(tpm))) stop("Check gene IDs and unique sample IDs.")
  gene_clean <- sub("\\..*", "", rownames(tpm))
  ord <- order(rowMeans(tpm), decreasing = TRUE)
  tpm <- tpm[ord, , drop = FALSE]
  gene_clean <- gene_clean[ord]
  keep <- !duplicated(gene_clean)
  tpm <- tpm[keep, , drop = FALSE]
  rownames(tpm) <- gene_clean[keep]

  gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
    keys = rownames(tpm), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  valid <- !is.na(gene_symbols)
  tpm_sym <- tpm[valid, , drop = FALSE]
  rownames(tpm_sym) <- gene_symbols[valid]
  tpm_sym <- tpm_sym[order(rowMeans(tpm_sym), decreasing = TRUE), , drop = FALSE]
  tpm_sym <- tpm_sym[!duplicated(rownames(tpm_sym)), , drop = FALSE]
  if (!nrow(tpm_sym)) stop("No Ensembl IDs mapped to gene symbols.")

  # 2. CIBERSORT on linear TPM, followed by PERK-UPR ssGSEA in primary tumors.
  tpm_sym_linear <- 2^tpm_sym - 1
  if (any(!is.finite(tpm_sym_linear))) stop("Check the log2(TPM + 1) input scale.")
  tumor_idx <- substr(colnames(tpm_sym), 14, 15) %in% "01"
  if (!any(tumor_idx)) stop("Expected TCGA primary-tumor sample barcodes (code 01).")
  tpm_tumor <- tpm_sym[, tumor_idx, drop = FALSE]
  clinical <- data.table::fread(args[2], sep = "\t", header = TRUE, data.table = FALSE)
  required <- c("sample", "OS", "OS.time")
  if (!all(required %in% names(clinical))) {
    stop("Survival input requires sample, OS, and OS.time columns.")
  }
  if (!any(colnames(tpm_tumor) %in% clinical$sample)) {
    stop("Expression and survival tables must use matching TCGA sample IDs.")
  }
  cibersort_res <- IOBR::deconvo_tme(eset = tpm_sym_linear,
    method = "cibersort", arrays = FALSE, perm = 1000)
  if (!all(c("ID", "Neutrophils_CIBERSORT") %in% names(cibersort_res))) {
    stop("Expected ID and Neutrophils_CIBERSORT columns from IOBR.")
  }
  if (anyDuplicated(cibersort_res$ID)) stop("CIBERSORT sample IDs must be unique.")

  gobp <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C5",
                          subcollection = "GO:BP")
  PERK_UPR_genes <- unique(gobp$gene_symbol[
    gobp$gs_name == "GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE"])
  found <- PERK_UPR_genes %in% rownames(tpm_tumor)
  if (!any(found)) stop("The PERK-UPR gene set must overlap the expression matrix.")
  message("PERK_UPR: matched ", sum(found), " of ", length(found), " genes.")
  ssgsea_param <- GSVA::ssgseaParam(exprData = tpm_tumor,
    geneSets = list(PERK_UPR = PERK_UPR_genes), normalize = TRUE)
  scores <- GSVA::gsva(ssgsea_param)
  neu_vec <- setNames(cibersort_res$Neutrophils_CIBERSORT, cibersort_res$ID)
  common_ids <- intersect(colnames(scores), names(neu_vec))
  if (!length(common_ids)) stop("CIBERSORT and ssGSEA sample IDs must match.")
  scores <- rbind(scores[, common_ids, drop = FALSE], Neu_sig = neu_vec[common_ids])
  score_df <- as.data.frame(t(scores))
  score_df$sample <- rownames(score_df)

  # 3. Original six-year follow-up and maximally selected cutpoints (minprop 0.3).
  surv_df <- prepare_survival_data(score_df, clinical)
  cp_neu <- survminer::surv_cutpoint(surv_df, time = "OS.time", event = "OS",
                                   variables = "Neu_sig", minprop = 0.3)
  surv_df$Neu_group <- factor(survminer::surv_categorize(cp_neu)$Neu_sig,
    levels = c("high", "low"), labels = c("High", "Low"))
  cp_perk <- survminer::surv_cutpoint(surv_df, time = "OS.time", event = "OS",
                                    variables = "PERK_UPR", minprop = 0.3)
  surv_df$PERK_group <- factor(survminer::surv_categorize(cp_perk)$PERK_UPR,
    levels = c("high", "low"), labels = c("High", "Low"))
  surv_df$combined_group <- factor(paste0(surv_df$PERK_group, "_", surv_df$Neu_group),
    levels = c("High_High", "High_Low", "Low_High", "Low_Low"))

  fit_neu <- survfit(Surv(OS.time, OS) ~ Neu_group, data = surv_df)
  p_neu <- survminer::ggsurvplot(fit_neu, data = surv_df,
    pval = TRUE, pval.method = TRUE, conf.int = FALSE, risk.table = TRUE,
    legend.title = "Neu", legend.labs = c("Neu_High", "Neu_Low"),
    palette = c("#E41A1C", "#4DAF4A"), xlab = "Time (months)",
    break.time.by = 12, xlim = c(0, 72), ggtheme = ggplot2::theme_bw())

  group_labels <- c(High_High = "PERK_High & Neu_High", High_Low = "PERK_High & Neu_Low",
                    Low_High = "PERK_Low & Neu_High", Low_Low = "PERK_Low & Neu_Low")
  group_colors <- c(High_High = "#E41A1C", High_Low = "#377EB8",
                    Low_High = "#4DAF4A", Low_Low = "#984EA3")
  present_groups <- levels(droplevels(surv_df$combined_group))
  fit_perk_neu <- survfit(Surv(OS.time, OS) ~ combined_group, data = surv_df)
  p_perk <- survminer::ggsurvplot(fit_perk_neu, data = surv_df,
    pval = TRUE, pval.method = TRUE, conf.int = FALSE, risk.table = TRUE,
    legend.title = "PERK & Neu", legend.labs = unname(group_labels[present_groups]),
    palette = unname(group_colors[present_groups]), xlab = "Time (months)",
    break.time.by = 12, xlim = c(0, 72), ggtheme = ggplot2::theme_bw())

  # 4. Export plots, scores, group assignments, cutpoints, and session details.
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (plot_name in c("KM_Neu", "KM_PERK_Neu_four_groups")) {
    grDevices::pdf(file.path(output_dir, paste0(plot_name, ".pdf")),
                   width = 10, height = 8, bg = "white")
    print(if (plot_name == "KM_Neu") p_neu else p_perk, newpage = FALSE)
    grDevices::dev.off()
  }
  write.csv(score_df, file.path(output_dir, "TCGA_LIHC_scores.csv"), row.names = FALSE)
  write.csv(surv_df, file.path(output_dir, "TCGA_LIHC_survival_groups.csv"), row.names = FALSE)
  cutpoints <- rbind(cp_neu$cutpoint, cp_perk$cutpoint)
  cutpoints$variable <- rownames(cutpoints)
  write.csv(cutpoints, file.path(output_dir, "cutpoints.csv"), row.names = FALSE)
  write.csv(data.frame(gene = PERK_UPR_genes, present = found),
            file.path(output_dir, "PERK_UPR_genes.csv"), row.names = FALSE)
  writeLines(c("seed=20261015", "CIBERSORT_permutations=1000", "minprop=0.3",
               "followup_months=72", "input_scale=log2(TPM+1)"),
             file.path(output_dir, "parameters.txt"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  message("Saved TCGA-LIHC survival analyses to ", output_dir)
}

if (sys.nframe() == 0L) main()
