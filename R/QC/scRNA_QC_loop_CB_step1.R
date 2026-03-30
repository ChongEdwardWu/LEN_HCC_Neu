### Section 0: Preparation --------------------------------------------------
# Clear the environment
rm(list = ls())
gc()

# Optional: load a personal R profile if you use one locally
if (file.exists("path_to_.radian_profile")) source("path_to_.radian_profile")

# Load required libraries
suppressPackageStartupMessages({
  library(rhdf5)                 # HDF5 I/O for CellBender .h5
  library(Matrix)                # Sparse matrices
  library(SingleCellExperiment)  # SCE container
  library(scuttle)               # QC helpers (addPerCellQC etc.)
  library(DropletUtils)          
  library(BiocParallel)
  library(scran)
  library(scater)
  library(ggplot2)
  library(gridExtra)
  library(scDblFinder)
  library(bluster)
  library(BiocSingular)
})

set.seed(123)
nworkers <- 16

### Section 1: Discover samples and loop through them ----------------------------------------
# 1. Specify the parent directory containing Cell Ranger count results
Basedir <- "path_to_project_root"
CBdir <- file.path(Basedir, "03_cellbender")

# 2. Get all subdirectory names (ignore non-directory files)
samples <- list.files(CBdir, pattern = "*", full.names = FALSE)
samples <- samples[sapply(file.path(CBdir, samples), dir.exists)]
# If you need to manually filter or sort, you can do so here, for example:
# samples <- samples[grepl("CTR|LEN", samples)]  # Example only

# 3. Loop through each subdirectory (i.e., each sample) to perform the analysis
for (sample in samples) {
  message("========== Analyzing sample: ", sample, " ==========")
  
  ### Section 2: Specify sample name and input -------------------------------------------
  # Create the working directory for this analysis (to store outputs and figures)
  # Change the base path as needed
  workdir <- file.path(Basedir, "04_R/QC", sample)
  if (!file.exists(workdir)) {
    dir.create(workdir, recursive = TRUE)
  }
  if (!file.exists(file.path(workdir, "figures"))) {
    dir.create(file.path(workdir, "figures"))
  }

  setwd(workdir)

  # Define the input path for the current sample's filtered_feature_bc_matrix
  inpath <- file.path(CBdir, sample, paste0(sample,"_filtered.h5"))
  message("Input cellbender path: ", inpath)

  # Read the 10X count matrix
  sce <- read10xCounts(inpath, type = "HDF5", col.names = TRUE)

  ## Section 3: Preprocess before QC -------------------------------------------
  # Deconvolution Normalization
  set.seed(123)
  clust <- quickCluster(sce, BPPARAM = MulticoreParam(nworkers))
  table(clust)
  sce <- computeSumFactors(sce,
    clusters = clust, min.mean = 0.1,
    BPPARAM = MulticoreParam(nworkers)
  )
  norm <- logNormCounts(sce)

  # Feature Selection
  dec <- modelGeneVar(norm)

  # Visualize the fit for feature selection
  fit <- metadata(dec)
  png(file = "figures/01_Feature_selection.png", width = 250, height = 250, units = "mm", res = 150)
  par(mfrow = c(1, 1))
  plot(fit$mean, fit$var,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression", main = ""
  )
  curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
  dev.off()

  # Select highly variable genes
  hvg.var <- getTopHVGs(dec, n = 2000)
  length(hvg.var)

  # Dimensionality Reduction and Clustering

  # PCA
  set.seed(123)
  norm <- fixedPCA(norm,
    subset.row = hvg.var,
    BSPARAM = RandomParam(), name = "PCA"
  )

  # UMAP
  norm <- runUMAP(norm,
    dimred = "PCA",
    BPPARAM = MulticoreParam(nworkers)
  )

  # Clustering
  colLabels(norm) <- clusterCells(norm,
    use.dimred = "PCA",
    BLUSPARAM = NNGraphParam(type = "jaccard")
  )

  # Save UMAP Clustering Plot
  sceumap <- plotReducedDim(norm, "UMAP", colour_by = "label", text_by = "label")
  ggsave(
    filename = "figures/02_Dimensionality_Reduction_and_Clustering.png",
    plot = sceumap,
    width = 250, height = 250, units = "mm",
    dpi = 150, device = "png", bg = "white"
  )

  # Doublet Detection
  set.seed(123)
  dbl <- scDblFinder(
    norm,
    clusters = colLabels(norm),
    nfeatures = 2000,
    BPPARAM = MulticoreParam(nworkers)
  )
  table(dbl$scDblFinder.class)
  norm$DoubletScore <- dbl$scDblFinder.score
  norm$DoubletClass <- dbl$scDblFinder.class

  ### Section 4: Preliminary QC ------------------------------------------------
  # Define gene markers
  is.mito <- grep("^(mt-|MT-)", rowData(norm)$Symbol)
  is.rp <- grep("^Rp[sl]|^RP[SL]", rowData(norm)$Symbol)
  is.hsp <- grep("^HSP|^DNAJ|^Hsp|^Dnaj", rowData(norm)$Symbol)
  is.hemoglobin <- grep("^Hb[ab]-|^HB[AB]", rowData(norm)$Symbol)
  is.plat <- grep("^Pf4|^Ppbp|^Itga2b|^Gp9|^Gp1ba|^Gp1bb|^Tubb1|^Vwf", rowData(norm)$Symbol)

  # Calculate QC metrics
  qc <- perCellQCMetrics(norm, subsets = list(
    Mito = is.mito, Heatshock = is.hsp, Rp = is.rp,
    RBC = is.hemoglobin, Plat = is.plat
  ))

  # Add QC metrics to colData
  colData(norm) <- cbind(colData(norm), qc)
  colnames(colData(norm))

  # (1) Filtering with fixed thresholds (modify thresholds as needed)
  qc.nexprs <- qc$detected < 200
  qc.mito <- qc$subsets_Mito_percent > 15
  qc.rbc <- qc$subsets_RBC_percent > 1
  qc.plat <- qc$subsets_Plat_percent > 1
  fix.discard <- qc.nexprs | qc.mito | qc.rbc | qc.plat
  table(fix.discard)
  norm$fix.discard <- fix.discard

  # Plot Mitochondrial % vs Detected Features
  png(file = "figures/03_Cell_Filtering.png", width = 250, height = 250, units = "mm", res = 150)
  plot(qc$detected, qc$subsets_Mito_percent,
    log = "x",
    xlab = "Detected features", ylab = "Mitochondrial %"
  )
  points(qc$detected[fix.discard], qc$subsets_Mito_percent[fix.discard],
    col = "dodgerblue", pch = 16, cex = 0.5
  )
  legend("topright",
    legend = c("Fix", "Adaptive"),
    col = c("dodgerblue", "red"), lty = 2, cex = 1
  )
  dev.off()

  # Define cells to discard based on fixed thresholds
  discard <- fix.discard
  norm$discard <- discard

  # Diagnostic Visualization
  diagnostic <- gridExtra::grid.arrange(
    plotColData(norm, x = "DoubletClass", y = "sum", colour_by = "DoubletClass") +
      ggtitle("Doublets"),
    plotColData(norm, x = "discard", y = "sum", colour_by = "discard") +
      ggtitle("Total count"),
    plotColData(norm,
      x = "discard", y = "subsets_Mito_percent",
      colour_by = "discard"
    ) + ggtitle("Mito percent"),
    plotColData(norm,
      x = "discard", y = "subsets_Rp_percent",
      colour_by = "discard"
    ) + ggtitle("Ribosomal protein percent"),
    plotColData(norm,
      x = "discard", y = "subsets_Heatshock_percent",
      colour_by = "discard"
    ) + ggtitle("Heatshock protein percent"),
    plotColData(norm, x = "sum", y = "detected", colour_by = "discard") +
      ggtitle("Detected features"),
    plotColData(norm,
      x = "sum", y = "subsets_RBC_percent",
      colour_by = "discard"
    ) + ggtitle("Hemoglobin percent"),
    plotColData(norm,
      x = "sum", y = "subsets_Plat_detected",
      colour_by = "discard"
    ) + ggtitle("Platelet percent"),
    ncol = 3
  )
  ggsave(
    filename = "figures/04_Diagnostic_Plots.png",
    plot = diagnostic,
    width = 500, height = 250, units = "mm",
    dpi = 150, device = "png", bg = "white"
  )

  # Cluster Quality Check
  clusQC <- gridExtra::grid.arrange(
    plotColData(norm, x = "label", y = "DoubletScore", colour_by = "DoubletClass") +
      ggtitle("Doublets"),
    plotColData(norm, x = "label", y = "sum", colour_by = "label") +
      ggtitle("Total count"),
    plotColData(norm, x = "label", y = "detected", colour_by = "label") +
      ggtitle("Detected features"),
    plotColData(norm, x = "label", y = "subsets_Mito_percent", colour_by = "label") +
      ggtitle("Mito percent"),
    plotColData(norm, x = "label", y = "subsets_Rp_percent", colour_by = "label") +
      ggtitle("Ribosomal genes percent"),
    plotColData(norm, x = "label", y = "subsets_Heatshock_percent", colour_by = "label") +
      ggtitle("Heatshock protein genes percent"),
    plotColData(norm, x = "label", y = "subsets_RBC_percent", colour_by = "label") +
      ggtitle("RBC percent"),
    plotColData(norm, x = "label", y = "subsets_Plat_percent", colour_by = "label") +
      ggtitle("Platelets percent"),
    plotColData(norm, x = "label", y = "discard", colour_by = "label") +
      ggtitle("Discard"),
    ncol = 3
  )
  ggsave(
    filename = "figures/05_Cluster_QC_Plots.png",
    plot = clusQC,
    width = 500, height = 500, units = "mm",
    dpi = 150, device = "png", bg = "white"
  )

  # View Filtering Results on UMAP
  sceumap <- gridExtra::grid.arrange(
    plotReducedDim(norm[, norm$discard == 0], "UMAP", colour_by = "label", text_by = "label") +
      ggtitle("Cells remained"),
    plotReducedDim(norm[, norm$discard == 1], "UMAP", colour_by = "label", text_by = "label") +
      ggtitle("Cells discarded"),
    ncol = 2
  )
  ggsave(
    filename = "figures/06_PreFiltering_UMAP.png",
    plot = sceumap,
    width = 250, height = 250, units = "mm",
    dpi = 150, device = "png", bg = "white"
  )

  ### Section 5: Save outputs --------------------------------------------------
  saveRDS(norm, file = paste0("01_QC_step1_", sample, ".rds"))

  message(paste0("Job ", sample, " has been successfully done at ", Sys.time(), "."))

  gc()
  graphics.off()

  # Proceed to the next sample
}
