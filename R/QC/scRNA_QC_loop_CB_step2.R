# -------- Section 0: Preparation ---------------------------------------------

# Clear workspace to avoid carry-over objects
rm(list = ls())
gc()

# Optional: load a personal R profile if you use one locally
if (file.exists("path_to_.radian_profile")) source("path_to_.radian_profile")
.libPaths()

suppressPackageStartupMessages({
    library(scater)
    library(edgeR)
    library(scran)
    library(BiocParallel)
    library(SingleR)
    library(scMCA)          # kept in case you re-enable later; not used below
    library(Seurat)
    library(SeuratDisk)
    library(SeuratWrappers)
    library(patchwork)
    library(ggplot2)
    library(stringr)
    library(tibble)
    library(Azimuth)
    library(matrixStats)
    library(Matrix)
    library(gridExtra)
    library(dplyr)
})
# httpgd::hgd()  # optional: interactive graphics device for plotting


## ---- Compatibility patch: Matrix ≥1.6 with older compiled binaries ---------
# Provide generic fallbacks to avoid method mismatches with newer Matrix.
if (!isGeneric("colSums")) setGeneric("colSums")
if (!isGeneric("rowMeans")) setGeneric("rowMeans")
if (!isGeneric("rowSums")) setGeneric("rowSums")

setMethod(
    "colSums", signature(x = "dgCMatrix"),
    function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...) {
        base::colSums(as.matrix(x), na.rm = na.rm)
    }
)
setMethod(
    "colSums", signature(x = "CsparseMatrix"),
    function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...) {
        base::colSums(as.matrix(x), na.rm = na.rm)
    }
)
setMethod(
    "rowMeans", signature(x = "dgCMatrix"),
    function(x, na.rm = FALSE, dims = 1L, ...) {
        base::rowMeans(as.matrix(x), na.rm = na.rm)
    }
)
setMethod(
    "rowSums", signature(x = "dgCMatrix"),
    function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...) {
        base::rowSums(as.matrix(x), na.rm = na.rm)
    }
)
setMethod(
    "rowSums", signature(x = "CsparseMatrix"),
    function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...) {
        base::rowSums(as.matrix(x), na.rm = na.rm)
    }
)
## ---------------------------------------------------------------------------


# -------- Section 1: User configuration ---------------------------------------

set.seed(123)
nworkers <- 20

## Species ("mm" for mouse, "hs" for human)
species <- "mm"

## Project layout (kept; personalized)
Basedir <- "path_to_project_root"

## Discover samples from QC step-1 output folders under 04_R/QC
QCbase <- file.path(Basedir, "04_R", "QC")
samples <- list.files(QCbase, full.names = FALSE)
samples <- samples[sapply(file.path(QCbase, samples), dir.exists)]

# Keep only folders that actually have step-1 RDS
has_step1 <- vapply(
    samples,
    \(s) file.exists(file.path(QCbase, s, paste0("01_QC_step1_", s, ".rds"))),
    logical(1)
)
samples <- samples[has_step1]

message(
    "Found ", length(samples), " sample(s) with QC step-1 output: ",
    paste(samples, collapse = ", ")
)


# -------- Main loop over samples ---------------------------------------------
for (sample in samples) {
    local({
        # Use local() so on.exit blocks are per-sample
        message("\n==== Processing sample: ", sample, " ====")

        # Paths that depend on 'sample' (kept)
        workdir <- file.path(QCbase, sample)
        dir.create(file.path(workdir, "figures"), recursive = TRUE, showWarnings = FALSE)
        setwd(workdir)

        step1_rds <- file.path(workdir, paste0("01_QC_step1_", sample, ".rds"))
        norm <- readRDS(step1_rds) # SingleCellExperiment with QC metrics

        # -------- Section 2: Identify unqualified clusters and finalize filtering ---

        # Mark clusters to discard (add numeric cluster labels if you have a list)
        colData(norm)$discard.cl <- ifelse(norm$label %in% c(
            # e.g., 5, 12  # put cluster labels to drop here if needed
        ), 1, 0)

        # Apply QC thresholds (BM defaults; adjust if needed)
        qc.nexprs  <- norm$detected < 100
        qc.mito    <- norm$subsets_Mito_percent > 10
        qc.rbc     <- norm$subsets_RBC_percent > 1
        qc.doublet <- norm$DoubletClass == "doublet"
        qc.cluster <- norm$discard.cl == 1

        # Final discard flag
        discard.new <- qc.nexprs | qc.mito | qc.rbc | qc.doublet | qc.cluster
        message("Summary of filtering:")
        print(table(discard.new))
        norm$discard <- discard.new

        # UMAP overview of kept vs. discarded cells
        sceumap <- gridExtra::grid.arrange(
            plotReducedDim(norm[, norm$discard == 0], "UMAP", colour_by = "label", text_by = "label") +
                ggtitle("Cells remained"),
            plotReducedDim(norm[, norm$discard == 1], "UMAP", colour_by = "label", text_by = "label") +
                ggtitle("Cells discarded"),
            ncol = 2
        )
        ggsave(
            filename = "figures/07_Final_Filtering_UMAP.png",
            plot = sceumap, width = 250, height = 250, units = "mm",
            dpi = 150, device = "png", bg = "white"
        )

        # Quick diagnostic: genes enriched in discarded pool
        rownames(rowData(norm)) <- rowData(norm)$Symbol
        lost      <- calculateAverage(counts(norm)[,  norm$discard])
        kept      <- calculateAverage(counts(norm)[, !norm$discard])
        logged    <- cpm(cbind(lost, kept), log = TRUE, prior.count = 2)
        logFC     <- logged[, 1] - logged[, 2]
        abundance <- rowMeans(logged)
        png(file = "figures/08_filtered_cells_DEGs.png", width = 250, height = 250, units = "mm", res = 150)
        plot(abundance, logFC, xlab = "Average count", ylab = "Log-FC (lost/kept)", pch = 16)
        dev.off()

        # Final kept SCE (post-filter)
        sce <- norm[, !norm$discard]
        sce$group <- sample

        # Unified CellID (keep raw barcodes in colnames; store unified ID in colData)
        barcodes    <- colnames(sce)
        sce$CellID  <- stringr::str_c(paste0(sample, ":"), stringr::str_sub(barcodes, 1, 16), "x")

        # ---- Normalize gene names to Gene Symbols for BOTH species ---------------
        # Fallback to gene ID if Symbol missing; enforce uniqueness
        sym   <- as.character(rowData(sce)$Symbol)
        gid   <- as.character(rowData(sce)$ID)
        use_id <- is.na(sym) | sym == ""
        sym[use_id] <- gid[use_id]
        sym <- make.unique(sym)
        rownames(sce) <- sym  # critical: ensures counts rownames are Symbols

        # Optional: mouse-only SingleR (ImmGen) annotation (kept; can be heavy)
        if (species == "mm") {
            load("path_to_immgen_reference_rdata")
            pred <- SingleR(
                test            = sce,
                ref             = immgen,
                labels          = immgen$label.fine,
                assay.type.test = 1,
                BPPARAM         = MulticoreParam(nworkers)
            )
            sce$CellType_immgen <- pred$pruned.labels
        }

        # Export SCE colData to RDS (for downstream joins in integration/figures)
        colData_tbl <- as_tibble(sce@colData)
        saveRDS(colData_tbl, file = paste0(sample, "_colData.rds"))


        # -------- Section 3: Build Seurat object from SCE RNA counts --------------
        # Always derive counts for Seurat from SCE (pure Gene Expression).

        cm <- counts(sce)                 # dgCMatrix in most pipelines
        rownames(cm) <- rownames(sce)     # already unified to Gene Symbol
        counts_mat <- as(cm, "dgCMatrix")
        source_tag <- "sce"

        # Safety checks
        stopifnot(nrow(counts_mat) > 2000, ncol(counts_mat) > 0)
        message("counts source = ", source_tag, "; dim = ", paste(dim(counts_mat), collapse = " x "))

        # Create Seurat object; genes are Symbols; cells are raw barcodes for now
        seu <- CreateSeuratObject(counts = counts_mat, project = sample)

        # Rename Seurat cells to the unified CellID and keep QC-passing cells only
        newnames <- stringr::str_c(paste0(sample, ":"), stringr::str_sub(Cells(seu), 1, 16), "x")
        seu <- RenameCells(seu, new.names = newnames)

        # Quick sanity check: are ALL Seurat cell IDs present in SCE colData?
        all_in <- all(Cells(seu) %in% colData_tbl$CellID)
        message("Checking CellID mapping: all(Cells(seu) %in% colData_tbl$CellID) = ", all_in)

        if (!all_in) {
            missing_ids <- setdiff(Cells(seu), colData_tbl$CellID)
            warning(sprintf(
                "CellIDs mismatch: %d Seurat cells not found in SCE colData for sample '%s'. Examples: %s",
                length(missing_ids), sample, paste(head(missing_ids, 10), collapse = ", ")
            ))
            # we're inside local({ ... }), so return to skip this sample safely
            return(invisible(NULL))
        }

        # (Optional) note if SCE has extra cells not present in Seurat
        extra_ids <- setdiff(colData_tbl$CellID, Cells(seu))
        if (length(extra_ids) > 0) {
            message(sprintf(
                "Note: %d CellIDs present in SCE colData but not in Seurat (they will be ignored).",
                length(extra_ids)
            ))
        }

        # Keep exactly the overlapping cells (order preserved)
        seu <- seu[, Cells(seu) %in% colData_tbl$CellID]


        # Merge SCE colData into Seurat meta.data (keeps QC metrics like mito%)
        if (!"CellID" %in% colnames(seu@meta.data)) {
            seu$CellID <- Cells(seu) # pre-create key column to avoid duplication
        }
        colData_tbl <- dplyr::distinct(colData_tbl, CellID, .keep_all = TRUE)

        meta_joined <- seu@meta.data %>%
            tibble::as_tibble() %>%
            dplyr::left_join(colData_tbl, by = "CellID")

        # Ensure row order matches Cells(seu) before assignment
        stopifnot(identical(meta_joined$CellID, Cells(seu)))
        seu@meta.data <- meta_joined %>% tibble::column_to_rownames("CellID")


        # -------- Section 4: Cell-cycle scoring and SCTransform -------------------

        # Use Seurat's updated CC gene sets; title-case for mouse
        cc <- Seurat::cc.genes.updated.2019
        to_mouse_case <- function(x) stringr::str_to_title(tolower(x))
        if (species == "mm") {
            s.features   <- intersect(to_mouse_case(cc$s.genes),   rownames(seu))
            g2m.features <- intersect(to_mouse_case(cc$g2m.genes), rownames(seu))
        } else {
            s.features   <- intersect(cc$s.genes,   rownames(seu))
            g2m.features <- intersect(cc$g2m.genes, rownames(seu))
        }

        # Normalize (RNA) → CellCycleScoring → compute CC.Diff for regression
        seu <- NormalizeData(seu, verbose = FALSE)
        seu <- CellCycleScoring(
            seu,
            s.features   = s.features,
            g2m.features = g2m.features,
            set.ident    = FALSE
        )
        seu$CCphase <- seu$Phase   # keep original Phase into CCphase
        seu$Phase   <- NULL
        seu$CC.Diff <- seu$S.Score - seu$G2M.Score

        # SCTransform with requested regressors (keep your configuration)
        vars_to_regress <- c("subsets_Mito_percent", "S.Score", "G2M.Score")
        seu <- SCTransform(
            seu,
            vst.flavor            = "v2",
            verbose               = TRUE,
            return.only.var.genes = FALSE,
            vars.to.regress       = vars_to_regress
        )
        # NOTE: No RunPCA here by design; integration will run PCA on shared features.


        # -------- Section 5: Optional Azimuth label transfer ----------------------
        # hs -> Human PBMC reference; mm -> Mouse PanSci reference (you provide files)

        # options(timeout = 10000)
        # az_ref_dir <- NULL
        # if (species == "hs") {
        #     az_ref_dir <- "path_to_azimuth_human_reference_dir"
        # } else if (species == "mm") {
        #     az_ref_dir <- "path_to_azimuth_mouse_reference_dir"
        # }

        # if (!is.null(az_ref_dir)) {
        #     ref_ok <- dir.exists(az_ref_dir) &&
        #               file.exists(file.path(az_ref_dir, "ref.Rds")) &&
        #               file.exists(file.path(az_ref_dir, "idx.annoy"))
        #     if (ref_ok) {
        #         message("Running Azimuth mapping using reference: ", az_ref_dir)
        #         az_res <- tryCatch(
        #             RunAzimuth(seu, reference = az_ref_dir),
        #             error = function(e) e
        #         )
        #         if (inherits(az_res, "error")) {
        #             warning("RunAzimuth failed for sample '", sample, "': ", conditionMessage(az_res))
        #         } else {
        #             # Append only new columns from Azimuth meta.data (robust handling)
        #             new_cols <- setdiff(colnames(az_res@meta.data), colnames(seu@meta.data))
        #             if (identical(rownames(seu@meta.data), rownames(az_res@meta.data))) {
        #                 if (length(new_cols)) {
        #                     seu@meta.data <- cbind(seu@meta.data,
        #                                            az_res@meta.data[, new_cols, drop = FALSE])
        #                 }
        #             } else {
        #                 md_old <- seu@meta.data; md_old$.__rn__ <- rownames(md_old)
        #                 md_new <- az_res@meta.data; md_new$.__rn__ <- rownames(md_new)
        #                 keep_new <- setdiff(colnames(md_new), colnames(md_old))
        #                 md_merge <- dplyr::left_join(md_old,
        #                                              md_new[, c(".__rn__", keep_new), drop = FALSE],
        #                                              by = ".__rn__")
        #                 rownames(md_merge) <- md_merge$.__rn__
        #                 md_merge$.__rn__ <- NULL
        #                 seu@meta.data <- md_merge
        #             }
        #             message("Azimuth columns added: ",
        #                     ifelse(length(new_cols), paste(new_cols, collapse = ", "), "none"))
        #         }
        #     } else {
        #         warning("Azimuth reference not found at ", az_ref_dir,
        #                 " (expecting files 'ref.Rds' and 'idx.annoy'); skipping Azimuth.")
        #     }
        # }


        # -------- Section 6: Save outputs -----------------------------------------

        saveRDS(seu, file = file.path(workdir, paste0("02_QC_step2_", sample, "_seu.rds")))
        message("✓ Done: ", sample, " at ", Sys.time(), " [source = ", source_tag, "]")

        gc()
    })
}
