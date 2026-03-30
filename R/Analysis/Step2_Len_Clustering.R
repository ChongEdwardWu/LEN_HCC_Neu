### ==========================================================================
### Step2: Clustering, annotation and marker/regulon analysis
### ==========================================================================

### ---------------- Section 0: Preparation / Global options ------------------

# WARNING: this removes all existing objects in the R session
rm(list = ls())
gc()

# Optional: load a personal R profile if you use one locally
if (file.exists("path_to_.radian_profile")) {
  suppressMessages(source("path_to_.radian_profile"))
}
# .libPaths()  # uncomment if you need to check library paths

# -----------------------------------------------------------------------------
# Working directory and output folders
# -----------------------------------------------------------------------------
workdir <- "path_to_data"
setwd(workdir)

if (!file.exists(file.path(workdir, "figures"))) {
  dir.create(file.path(workdir, "figures"))
}
if (!file.exists(file.path(workdir, "results"))) {
  dir.create(file.path(workdir, "results"))
}

# -----------------------------------------------------------------------------
# Packages
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(SeuratWrappers)
  library(SeuratObject)
  library(patchwork)
  library(ggplot2)
  library(future)
  library(SingleR)
  library(BiocParallel)
  library(openxlsx2)
  # Optional analysis blocks:
  # library(GSEABase)
  # library(AUCell)
  # library(clusterProfiler)
  # library(org.Mm.eg.db)
})

# httpgd::hgd()  # enable httpgd if you want interactive plotting in browser

# -----------------------------------------------------------------------------
# Parallelization and global options
# -----------------------------------------------------------------------------
set.seed(123)

nworkers <- 20

# First set a safe default (sequential), then override with multicore
plan("sequential")
plan("multicore", workers = nworkers)

# Species flag (used e.g. for blacklist; "mm" for mouse, "hs" for human)
species <- "mm"   # choices: "mm" or "hs"

# -----------------------------------------------------------------------------
# Input: integrated Seurat object with pySCENIC assays
# -----------------------------------------------------------------------------
integrated_rds <- "results/01.5_LEN_pySCENIC_results.rds"
stopifnot(file.exists(integrated_rds))

seu <- readRDS(integrated_rds)
message("Loaded integrated object: cells = ", ncol(seu),
        "  genes = ", nrow(seu))

### ---------------- Section 1: Clustering handle (@res=1) --------------------
# Here we *reuse* the clustering done during integration.
# You only need to choose which resolution to work with.

res <- 1  # resolution you decided to use in previous step

# UMAP split by group (treatment / condition) colored by cluster
DimPlot(
  seu,
  reduction = "umap",
  split.by  = "group",
  group.by  = paste0("integrated_snn_res.", res),
  label     = TRUE
) + coord_fixed(ratio = 1)

# Set the active identity and a coarse level-1 CellType factor
Idents(seu) <- seu$CellType_l1 <- factor(
  seu@meta.data[[paste0("integrated_snn_res.", res)]]
)

# If CCphase is missing but Seurat's Phase is present, reuse it
if (!"CCphase" %in% colnames(seu@meta.data) && "Phase" %in% colnames(seu@meta.data)) {
  seu$CCphase <- seu$Phase
}

### ---------------- Section 2: Automated annotations -------------------------
# Here we gather automated annotations (SingleR, AUCell, MCA, PanSci, etc.)
# into a single meta column 'annot_auto' for quick reference.
# SingleR block is kept as template but commented out.

DefaultAssay(seu) <- "RNA"

# -----------------------------------------------------------------------------
# 2.1 SingleR + ImmGen (optional; robust for mouse hematopoiesis)
# -----------------------------------------------------------------------------
# immgen_rdata <- "path_to_immgen_reference_rdata"
# if (file.exists(immgen_rdata)) {
#   load(immgen_rdata)  # loads 'immgen'
#   message("ImmGen reference loaded.")
# 
#   DefaultAssay(seu) <- "RNA"
#   seu <- NormalizeData(seu, verbose = FALSE)
#   sce <- as.SingleCellExperiment(seu)
# 
#   pred <- SingleR(
#     test            = sce,
#     ref             = immgen,
#     labels          = immgen$label.fine,
#     assay.type.test = 1,
#     BPPARAM         = MulticoreParam(nworkers)
#   )
# 
#   # Fine-grained ImmGen labels for each cell
#   seu$CellType_immgen <- pred$pruned.labels
#   # Optional diagnostic plot:
#   # plotScoreHeatmap(pred)
# } else {
#   warning("ImmGen RData not found: ", immgen_rdata,
#           " — skipping SingleR ImmGen annotation.")
# }

# -----------------------------------------------------------------------------
# 2.2 AUCell (optional example)
# -----------------------------------------------------------------------------
# Example template for building AUCell-based signatures:
# ClMk <- read.csv(file = "HSPC cluster marker genes.csv")
# MkClusters <- unique(ClMk$Cluster.annotation)
# ClMarkers <- vector("list", length(MkClusters))
# names(ClMarkers) <- MkClusters
# for (i in seq_along(MkClusters)) {
#   this_cluster <- MkClusters[i]
#   ClMarkers[[this_cluster]] <- as.character(
#     ClMk$Gene.symbol[ClMk$Cluster.annotation %in% this_cluster]
#   )
# }
# save(ClMarkers, file = "HSPC_cluster_marker_genes.RData")

# -----------------------------------------------------------------------------
# 2.3 Combine available auto-annotations into a single column
# -----------------------------------------------------------------------------
auto_cols <- intersect(
  c("CellType_immgen", "annot_HSPC_AUCell", "CellType_mca", "CellType_pansci"),
  colnames(seu@meta.data)
)

if (length(auto_cols) > 0) {
  seu$annot_auto <- factor(
    apply(
      seu@meta.data[, auto_cols, drop = FALSE],
      1,
      function(v) paste(na.omit(as.character(v)), collapse = "; ")
    )
  )
}

# -----------------------------------------------------------------------------
# 2.4 Summarize cluster annotation composition
# -----------------------------------------------------------------------------
cluster_annot <- tibble(
  Cluster  = seu$seurat_clusters,
  Source   = seu$group,
  CellType = seu$annot_auto
) %>%
  group_by(Cluster, CellType) %>%
  summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(Cluster) %>%
  mutate(
    total.no = sum(no.cell),
    perc     = 100 * no.cell / total.no
  ) %>%
  arrange(Cluster, dplyr::desc(perc)) %>%
  dplyr::slice_max(order_by = perc, n = 5)

# View(cluster_annot)  # interactive inspection if needed

# -----------------------------------------------------------------------------
# 2.5 Cluster composition across groups
# -----------------------------------------------------------------------------
source_cluster <- tibble(
  Cluster = seu$seurat_clusters,
  Group   = seu$group
) %>%
  group_by(Cluster, Group) %>%
  summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(Group) %>%
  mutate(
    total.no = sum(no.cell),
    perc     = 100 * no.cell / total.no
  ) %>%
  dplyr::select(Cluster, Group, perc)

# Quick barplot: cluster distribution by group
ggplot(source_cluster, aes(x = Group, y = perc, fill = Cluster)) +
  geom_col(colour = "black") +
  coord_fixed(ratio = 1 / 10) +
  theme_bw() +
  xlab("Group") +
  ylab("Cell fraction (%)")

# -----------------------------------------------------------------------------
# 2.6 Cell cycle phase composition per cluster
# -----------------------------------------------------------------------------
source_CCphase <- tibble(
  Cluster = seu$seurat_clusters,
  CCphase = seu$CCphase
) %>%
  group_by(Cluster, CCphase) %>%
  summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(CCphase) %>%
  mutate(
    total.no = sum(no.cell),
    perc     = 100 * no.cell / total.no
  ) %>%
  dplyr::select(Cluster, CCphase, perc)

# -----------------------------------------------------------------------------
# 2.7 Export summary tables to Excel
# -----------------------------------------------------------------------------
save_file_name <- "results/Step2_Len_cluster_celltype.xlsx"
wb <- wb_workbook()

# Sheet 1: cluster_annot
wb <- wb_add_worksheet(wb, sheet = "cluster_annot")
wb <- wb_add_data(
  wb,
  sheet = "cluster_annot",
  x     = as.data.frame(cluster_annot)
)

# Sheet 2: source_cluster (wide format)
wb <- wb_add_worksheet(wb, sheet = "source_cluster")
wb <- wb_add_data(
  wb,
  sheet = "source_cluster",
  x     = as.data.frame(source_cluster)
)

# Sheet 3: source_CCphase (wide format)
source_CCphase_wide <- source_CCphase %>%
  tidyr::pivot_wider(names_from = Cluster, values_from = perc)

wb <- wb_add_worksheet(wb, sheet = "source_CCphase")
wb <- wb_add_data(
  wb,
  sheet = "source_CCphase",
  x     = as.data.frame(source_CCphase_wide)
)

wb_save(wb, save_file_name, overwrite = TRUE)

### ---------------- Section 3: Lineage marker dotplot (CellType_l1 helper) ---

Idents(seu) <- "seurat_clusters"

# Compact mouse hematopoiesis / immune lineage marker panel
# Reference: Kucinski et al., Cell Stem Cell 2024 (adapted)
lineage_panel_mm <- list(
  HSC                         = c("Procr", "Mecom"),
  `Megakaryocyte (progenitor)` = c("Pf4", "Vwf", "Itga2b", "Gp9", "Itgb3"),
  `Erythroid (progenitor)`     = c("Klf1", "Gata1", "Hba-a1"),
  `Neu (progenitor)`           = c("Elane", "Prtn3", "Ms4a3"),
  `Basophil (progenitor)`      = c("Prss34", "Mcpt8"),
  `Eosinophil (progenitor)`    = c("Prg2", "Prg3", "Epx"),
  `Mast cell (progenitor)`     = c("Kit", "Cma1", "Gzmb"),
  `Mono/DC progenitor`         = c("Csf1r", "Mpo"),
  `Lymphoid progenitor`        = c("Il7r", "Dntt"),
  `T cell (progenitor)`        = c("Cd3e", "Cd8a", "Bcl11b"),
  `B cell (progenitor)`        = c("Ly6d", "Vpreb3", "Cd79a"),
  ILC                          = c("Gata3", "Id2", "Il7r"),
  Mono                         = c("Cd14", "Ctsg", "Ly6c1", "Ly6c2"),
  pDC                          = c("Irf8", "Siglech", "Ly6d"),
  `Ifn-activated cells`        = c("Ifitm3", "Irf7", "Isg15"),
  `Complement expressing`      = c("C1qa", "C1qb"),
  DCs                          = c("H2-Aa", "Xcr1", "Itgax")
)

# Keep only genes that are present in the object
gene_panel <- unique(unlist(lineage_panel_mm))
gene_panel <- intersect(gene_panel, rownames(seu))

# DotPlot of lineage markers by seurat clusters (SCT assay)
DotPlot(
  seu,
  group.by = "seurat_clusters",
  assay    = "SCT",
  features = gene_panel
) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.title       = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    legend.direction = "vertical",
    legend.position  = "bottom"
  ) +
  scale_color_gradientn(
    colours = c("white", "#fde3d8", "#ed684e", "#d21e20", "#981917")
  )

### ---------------- Section 4: Manual lineage annotation ----------------------

# Manual mapping from Seurat clusters to broad lineages (CellType_l1).
# Update these labels if cluster identities change.
seu$CellType_l1 <- factor(dplyr::recode(
  as.character(seu$seurat_clusters),
  "1"  = "B",
  "2"  = "Neu",
  "3"  = "Neu",
  "4"  = "Neu",
  "5"  = "T",
  "6"  = "T",
  "7"  = "Mono",
  "8"  = "Neu",
  "9"  = "T",
  "10" = "NK",
  "11" = "T",
  "12" = "Mac",
  "13" = "T",
  "14" = "DC",
  "15" = "T",
  "16" = "pDC",
  "17" = "Baso",
  "18" = "B"
), levels = c("T", "B", "NK", "DC", "pDC", "Mono", "Mac", "Neu", "Baso"))

# Fine-grained manual mapping used by downstream neutrophil analyses.
seu$CellType_l2 <- dplyr::recode(
  as.character(seu$seurat_clusters),
  "5"  = "c1_T_CD8eff_Cxcr3",
  "6"  = "c2_T_CD8exh_Pdcd1Lag3",
  "9"  = "c3_T_Treg_Foxp3_Il2ra",
  "11" = "c4_T_Naive_Tcf7_Il7r",
  "13" = "c5_T_CD8cycle_Mki67",
  "15" = "c6_T_NKT_Zbtb16_Tcf7",
  "1"  = "c7_B_Fo_Cd79a",
  "18" = "c8_B_Plasma_Jchain",
  "10" = "c9_NK_Eomes_Gzmb",
  "14" = "c10_DC_cDC1_Xcr1",
  "16" = "c11_pDC_Siglech_Tcf4",
  "7"  = "c12_Mono_Ly6C_Csf1r",
  "12" = "c13_Mac_C1qa_Mrc1",
  "2"  = "c14_Neu_ISG_Isg15",
  "3"  = "c15_Neu_S100a8a9",
  "4"  = "c16_Neu_Spp1_Ccl3",
  "8"  = "c17_Neu_Mature_Camp",
  "17" = "c18_Baso_Mcpt8_Il4"
)

seu$CellType_l2 <- factor(
  seu$CellType_l2,
  levels = c(
    "c1_T_CD8eff_Cxcr3",
    "c2_T_CD8exh_Pdcd1Lag3",
    "c3_T_Treg_Foxp3_Il2ra",
    "c4_T_Naive_Tcf7_Il7r",
    "c5_T_CD8cycle_Mki67",
    "c6_T_NKT_Zbtb16_Tcf7",
    "c7_B_Fo_Cd79a",
    "c8_B_Plasma_Jchain",
    "c9_NK_Eomes_Gzmb",
    "c10_DC_cDC1_Xcr1",
    "c11_pDC_Siglech_Tcf4",
    "c12_Mono_Ly6C_Csf1r",
    "c13_Mac_C1qa_Mrc1",
    "c14_Neu_ISG_Isg15",
    "c15_Neu_S100a8a9",
    "c16_Neu_Spp1_Ccl3",
    "c17_Neu_Mature_Camp",
    "c18_Baso_Mcpt8_Il4"
  )
)

table(seu$seurat_clusters, seu$CellType_l1)
table(seu$group, seu$CellType_l1)
table(seu$seurat_clusters, seu$CellType_l2)

DimPlot(
  seu,
  reduction = "umap",
  group.by  = "CellType_l1",
  label     = TRUE
) + coord_fixed(1)

DimPlot(
  seu,
  reduction = "umap",
  group.by  = "CellType_l1",
  split.by  = "group",
  label     = TRUE
) + coord_fixed(1)

### ---------------- Section 5: Conserved markers and regulon markers ---------
save_file_name <- "results/Step2_Len_cluster_annot_final.xlsx"
wb <- wb_workbook()

# Sheet 1: cluster_annot
cluster_annot <- tibble(
  Cluster  = seu$CellType_l2,
  Source   = seu$group,
  CellType = seu$annot_auto
) %>%
  group_by(Cluster, CellType) %>%
  summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(Cluster) %>%
  mutate(
    total.no = sum(no.cell),
    perc     = 100 * no.cell / total.no
  ) %>%
  arrange(Cluster, dplyr::desc(perc)) %>%
  dplyr::slice_max(order_by = perc, n = 5)
wb <- wb_add_worksheet(wb, sheet = "cluster_annot")
wb <- wb_add_data(
  wb,
  sheet = "cluster_annot",
  x     = as.data.frame(cluster_annot)
)

# Sheet 2: source_cluster
source_cluster <- tibble(
  Cluster = seu$CellType_l2,
  Group   = seu$group
) %>%
  group_by(Cluster, Group) %>%
  summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(Group) %>%
  mutate(
    total.no = sum(no.cell),
    perc     = 100 * no.cell / total.no
  )

wb <- wb_add_worksheet(wb, sheet = "source_cluster")
wb <- wb_add_data(
  wb,
  sheet = "source_cluster",
  x     = as.data.frame(source_cluster)
)

# Build an RNA feature whitelist for conserved marker testing.
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)
rna_genes <- rownames(seu[["RNA"]])

build_bad_features <- function(gene_names) {
  hist_genes <- grep("^Hist", gene_names, ignore.case = TRUE, value = TRUE)
  hb_genes <- grep("^Hb[ab]-|^HB(?!P)", gene_names, perl = TRUE, value = TRUE)
  mt_genes <- grep("^mt-|^MT-", gene_names, ignore.case = TRUE, value = TRUE)
  rps_genes <- grep("^Rp[sl]|^RP[SL]", gene_names, ignore.case = TRUE, value = TRUE)
  rik_genes <- grep("Rik$", gene_names, value = TRUE)
  pseudo_genes <- grep("-(ps|ps[0-9]+)$", gene_names, ignore.case = TRUE, value = TRUE)
  gencode_genes <- grep("^Gm[0-9]+", gene_names, value = TRUE)
  unique(c(hist_genes, hb_genes, mt_genes, rps_genes, rik_genes, pseudo_genes, gencode_genes))
}

features_rna <- setdiff(rna_genes, build_bad_features(rna_genes))
conserved_df <- NULL
if ("group" %in% colnames(seu@meta.data) && length(unique(seu$group)) > 1) {
  Idents(seu) <- "CellType_l2"
  cluster_ids <- levels(seu$CellType_l2)
  conserved_list <- lapply(cluster_ids, function(cl) {
    tryCatch({
      FindConservedMarkers(
        seu,
        ident.1              = cl,
        grouping.var         = "group",
        only.pos             = TRUE,
        min.pct              = 0.3,
        min.cells.per.ident  = 10,
        logfc.threshold      = 0.25,
        test.use             = "wilcox",
        features             = features_rna
      ) %>% dplyr::mutate(cluster = cl, gene = rownames(.))
    }, error = function(e) {
      warning("Conserved markers failed for cluster ", cl, ": ", e$message)
      NULL
    })
  })
  conserved_df <- conserved_list %>% bind_rows() %>%
    group_by(cluster) %>% arrange(max_pval, .by_group = TRUE)

  wb <- wb_add_worksheet(wb, "Conserved_all")
  wb <- wb_add_data(wb, "Conserved_all", as.data.frame(conserved_df))
  wb_save(wb, save_file_name, overwrite = TRUE)
} else {
  message("Skip conserved markers: 'group' not found or only one level present.")
}


# -----------------------------------------------------------------------------
# 5.5 Regulon markers by cluster (AUC assay)
# -----------------------------------------------------------------------------
Idents(seu) <- "CellType_l2"
seu_AUC <- FindAllMarkers(
  seu,
  assay           = "AUC",
  only.pos        = TRUE,
  min.pct         = 0.1,
  logfc.threshold = 0.25
)

seu_AUC <- seu_AUC %>%
  group_by(cluster) %>%
  arrange(p_val_adj, .by_group = TRUE)

wb <- wb_add_worksheet(wb, "Clus_Regulon_AUC")
wb <- wb_add_data(wb, "Clus_Regulon_AUC", as.data.frame(seu_AUC))
wb_save(wb, save_file_name, overwrite = TRUE)

# -----------------------------------------------------------------------------
# 5.6 Regulon markers by cluster (Bin assay: on/off)
# -----------------------------------------------------------------------------
seu_Bin <- FindAllMarkers(
  seu,
  assay           = "Bin",
  only.pos        = TRUE,
  min.pct         = 0.1,
  logfc.threshold = 0.25
)

seu_Bin <- seu_Bin %>%
  group_by(cluster) %>%
  arrange(p_val_adj, .by_group = TRUE)

wb <- wb_add_worksheet(wb, "Clus_Regulon_Bin")
wb <- wb_add_data(wb, "Clus_Regulon_Bin", as.data.frame(seu_Bin))
wb_save(wb, save_file_name, overwrite = TRUE)

### ---------------- Section 6: Save outputs ----------------------------------
saveRDS(seu, file = "results/02_Annotation_LSK.rds")
writeLines(capture.output(sessionInfo()), "session_info.txt")
