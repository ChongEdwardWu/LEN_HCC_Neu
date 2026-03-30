### Step 1: Data integration
### Seurat SCT-RPCA integration, dimension reduction, and clustering


### Section 0: Preparation --------------------------------------------------
rm(list = ls()); gc()

# Optional: load a personal R profile if you use one locally
if (file.exists("path_to_.radian_profile")) source("path_to_.radian_profile")
# .libPaths()  # keep your personal lib setup if needed

# Set your working directory
workdir <- "path_to_data"
setwd(workdir)

# Create directories for figures/results (idempotent)
if (!file.exists(file.path(workdir, "figures"))) dir.create(file.path(workdir, "figures"))
if (!file.exists(file.path(workdir, "results"))) dir.create(file.path(workdir, "results"))

### Section 1: Data integration ---------------------------------------------

## Using Seurat to integrate datasets ##
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(patchwork)
  library(tidyverse)
  library(ggplot2)
  library(clustree)
})

set.seed(123)
nworkers <- 8
# httpgd::hgd()

## Species configuration
species <- "mm"  # "mm" for mouse, "hs" for human

## Discover samples that have QC Step-2 outputs
QCdir   <- file.path(workdir, "QC")
samples <- list.files(QCdir, full.names = FALSE)
samples <- samples[sapply(file.path(QCdir, samples), dir.exists)]

# Keep only folders that actually have the *Step-2* RDS (updated name)
has_step2 <- vapply(
  samples,
  \(s) file.exists(file.path(QCdir, s, paste0("02_QC_step2_", s, "_seu.rds"))),
  logical(1)
)
samples <- samples[has_step2]

message("Integration samples (QC Step-2 detected): ", paste(samples, collapse = ", "))

# ---- Load Step-2 Seurat objects into a named list --------------------------
Seu_list <- list()
for (s in samples) {
  fp <- file.path(QCdir, s, paste0("02_QC_step2_", s, "_seu.rds"))
  obj <- readRDS(fp)
  # Optional sanity checks (all SCT; metadata carried through)
  stopifnot(DefaultAssay(obj) %in% c("SCT", "RNA"))  # SCT expected; RNA also fine
  Seu_list[[s]] <- obj
  rm(obj)
}
names(Seu_list) <- samples

# ---- Feature selection (strict common set for RPCA) ------------------------
features <- SelectIntegrationFeatures(object.list = Seu_list, nfeatures = 2000)

# Important: use only genes present in *all* objects to avoid zero-length combines
present_in_all <- Reduce(intersect, lapply(Seu_list, rownames))
features <- intersect(features, present_in_all)
message("Integration features after intersection: ", length(features))

# Precompute Pearson residuals where needed
Seu_list <- PrepSCTIntegration(
  object.list     = Seu_list,
  anchor.features = features,
  verbose         = TRUE
)

# RPCA requires PCA computed on the *same* features
Seu_list <- lapply(
  Seu_list,
  function(x) RunPCA(x, features = features, npcs = 50, verbose = FALSE)
)

# ---- Anchors & Integration (SCT-RPCA) --------------------------------------
# Choose explicit references if present; else let Seurat decide
seurat_ref <- c("CTR","LEN")
# seurat_ref <- c()
ref_idx    <- match(seurat_ref, names(Seu_list))
ref_idx    <- ref_idx[!is.na(ref_idx)]
if (!length(ref_idx)) ref_idx <- NULL

anchors <- FindIntegrationAnchors(
  object.list          = Seu_list,
  normalization.method = "SCT",
  reduction            = "rpca",
  anchor.features      = features,
  dims                 = 1:30,
  reference            = ref_idx,
  verbose              = TRUE
)

seu <- IntegrateData(
  anchorset            = anchors,
  normalization.method = "SCT",
  dims                 = 1:30,
  verbose              = TRUE
)

# Save memory
rm(Seu_list); gc()


### Section 2: Dimension reduction and clustering ----------------------------
## Standard workflow for visualization and clustering on the integrated assay

# Exclude unwanted genes from PCA (cell cycle, ribosomal, etc.)
# Cell-cycle genes (human list -> title-case for mouse)
cc_genes_all <- unique(unlist(Seurat::cc.genes.updated.2019))
to_mouse_case <- function(x) stringr::str_to_title(tolower(x))
if (species == "mm") {
  cc_genes <- to_mouse_case(cc_genes_all)
} else {
  cc_genes <- cc_genes_all
}

rna_genes <- rownames(seu[["RNA"]])
# Histone
hist_genes   <- grep("^Hist", rna_genes, ignore.case = TRUE, value = TRUE)
# Hemoglobin
hb_genes     <- grep("^Hb[ab]-|^HB(?!P)", rna_genes, perl = TRUE, value = TRUE)
# Mitochondrial
mt_genes     <- grep("^mt-|^MT-", rna_genes, ignore.case = TRUE, value = TRUE)
# Ribosomal
rps_genes    <- grep("^Rp[sl]|^RP[SL]", rna_genes, ignore.case = TRUE, value = TRUE)
# Riken (mouse)
rik_genes     <- grep("Rik$", rna_genes, value = TRUE)
# Alu elements
if (species == "hs") {
  alu_genes <- grep("^AL[0-9]+", rna_genes, value = TRUE)
} else {
  alu_genes <- character(0)
}
# Pseudogenes
pseudo_genes  <- grep("-(ps|ps[0-9]+)$", rna_genes, ignore.case = TRUE, value = TRUE)
# microRNAs
mir_genes    <- grep("^Mir", rna_genes, ignore.case = TRUE, value = TRUE)
# Gencode-style “Gm” (mouse models)
gencode_genes <- grep("^Gm[0-9]+", rna_genes, value = TRUE)

# Combine
bad_features <- unique(c(
  # cc_genes,
  hist_genes, hb_genes, rps_genes, mt_genes,
  rik_genes, gencode_genes, pseudo_genes, alu_genes
))

# Work on the integrated assay from here
DefaultAssay(seu) <- "integrated"

# PCA features: integrated HVGs minus “bad” features
PCA_features <- setdiff(seu[["integrated"]]@var.features, bad_features)
length(PCA_features)

seu <- RunPCA(
  seu,
  npcs     = 50,
  verbose  = TRUE,
  features = PCA_features
)

# Elbow/JackStraw (optional)
ElbowPlot(seu, ndims = 50)
# seu <- JackStraw(seu, num.replicate = 100)
# seu <- ScoreJackStraw(seu, dims = 1:20)
# JackStrawPlot(seu, dims = 1:20)

## UMAP/Neighbors
suppressPackageStartupMessages(library(future))
plan("sequential")
# plan("multicore", workers = nworkers)  # uncomment if you prefer
set.seed(123)
DefaultAssay(seu) <- "integrated"

# 
pca_dims_use <- 1:20
seu <- FindNeighbors(
  seu,
  reduction = "pca", dims = pca_dims_use,
  k.param = 50, prune.SNN = 0
)
seu <- RunUMAP(
  seu,
  reduction = "pca", dims = pca_dims_use,
  n.neighbors = 80, min.dist = 0.6, spread = 1.3,
  metric = "cosine", seed.use = 123
)


## Estimate clustering across multiple resolutions
# seu@meta.data[,grep("integrated_snn_res.",colnames(seu@meta.data),value = TRUE)] <- NULL
# seu <- FindClusters(seu, resolution = 1)  # single-res example

seu <- FindClusters(
  seu,
  resolution = c(seq(0.1, 1, 0.1))
)

# Relevel cluster names (avoid “0” label). Keep this pattern you use:
for (i in seq_along(grep("integrated_snn_res.", colnames(seu@meta.data)))) {
  j <- grep("integrated_snn_res.", colnames(seu@meta.data))[i]
  k <- seu@meta.data[, j]
  levels(k) <- as.character(seq_len(length(levels(k))))
  seu@meta.data[, j] <- k
}

# Visualize cluster hierarchy across resolutions
suppressPackageStartupMessages(library(clustree))
clustree(seu@meta.data, prefix = "integrated_snn_res.", return = "plot")

# ! resolution choice
res <- 1
# Plot 1: group difference in clusters (UMAP split)
DimPlot(
  seu,
  reduction = "umap",  #split.by = "group",
  group.by = paste0("integrated_snn_res.", res), label = TRUE
) + coord_fixed(ratio = 1)

# ! set resolution
Idents(seu) <- seu$seurat_clusters <- factor(seu@meta.data[[paste0("integrated_snn_res.", res)]])
# seu@meta.data[,grep("integrated_snn_res.",colnames(seu@meta.data),value = TRUE)] <- NULL

# QC overview (these metrics were merged from SCE during Step-2)
VlnPlot(
  seu,
  features = c("sum", "detected", "subsets_Mito_percent", "subsets_Rp_percent", "subsets_Heatshock_percent"),
  pt.size = 0
)

# ---------------------------------------------------------------------------
# Cluster summaries and simple exports
# NOTE: Step-2 moved cell-cycle scoring to Seurat (Phase column). If CCphase
# is absent (from scran::cyclone), fall back to Phase to keep your plots working.
if (!"CCphase" %in% colnames(seu@meta.data) && "Phase" %in% colnames(seu@meta.data)) {
  seu$CCphase <- seu$Phase
}

# View cluster annotation (ImmGen labels should exist if Step-2 added them)
cluster_annot <- tibble(
  seu$seurat_clusters,
  seu$group,
  seu$CellType_immgen
) %>%
  set_names("Cluster", "Source", "CellType") %>%
  group_by(Cluster, CellType) %>%
  summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(Cluster) %>%
  mutate(
    total.no = sum(no.cell),
    perc     = 100 * no.cell / total.no
  ) %>%
  arrange(Cluster, dplyr::desc(perc)) %>%
  dplyr::slice_head(n = 5)

# View and plot cluster distribution
source_cluster <- tibble(
  seu$seurat_clusters,
  seu$group
) %>%
  set_names("Cluster", "Group") %>%
  group_by(Cluster, Group) %>%
  summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(Group) %>%
  mutate(
    total.no = sum(no.cell),
    perc     = 100 * no.cell / total.no
  ) %>%
  dplyr::select(Cluster, Group, perc)

ggplot(source_cluster, aes(x = Group, y = perc, fill = Cluster)) +
  geom_col(colour = "black") +
  # scale_fill_manual(values = cl) +  # keep for custom palettes
  coord_fixed(ratio = 1 / 10) +
  theme_bw() +
  xlab("group") +
  ylab("%")

# View and plot CCphase distribution (uses fallback to Phase if needed)
source_CCphase <- tibble(
  seu$seurat_clusters,
  seu$CCphase
) %>%
  set_names("Cluster", "CCphase") %>%
  group_by(Cluster, CCphase) %>%
  summarise(no.cell = n(), .groups = "drop_last") %>%
  group_by(CCphase) %>%
  mutate(
    total.no = sum(no.cell),
    perc     = 100 * no.cell / total.no
  ) %>%
  dplyr::select(Cluster, CCphase, perc)

# Export results (xlsx)
suppressPackageStartupMessages(library(openxlsx2))

save_file_name <- "results/Step1_cluster_celltype.xlsx"
wb <- wb_workbook()

wb <- wb_add_worksheet(wb, sheet = "cluster_annot")
wb <- wb_add_data(wb, sheet = "cluster_annot", x = as.data.frame(cluster_annot))

source_cluster_wide <- source_cluster %>%
  tidyr::pivot_wider(names_from = Cluster, values_from = perc)
wb <- wb_add_worksheet(wb, sheet = "source_cluster")
wb <- wb_add_data(wb, sheet = "source_cluster", x = as.data.frame(source_cluster_wide))

source_CCphase_wide <- source_CCphase %>%
  tidyr::pivot_wider(names_from = Cluster, values_from = perc)
wb <- wb_add_worksheet(wb, sheet = "source_CCphase")
wb <- wb_add_data(wb, sheet = "source_CCphase", x = as.data.frame(source_CCphase_wide))

wb_save(wb, save_file_name, overwrite = TRUE)

### Section 3: Filter bad cells/clusters (optional) ---------------------------
# ## Example sketch (keep these commented lines for quick manual filtering):
# discardCl <- seu$seurat_clusters %in% c(18)
# discardDim <- (((seu[["umap"]]@cell.embeddings[,1] > 3 | seu[["umap"]]@cell.embeddings[,1] < -2) | (seu[["umap"]]@cell.embeddings[,2] > 7 | seu[["umap"]]@cell.embeddings[,2] < 0)) & (seu$seurat_clusters %in% c(8)))
# discardMito <- seu$subsets_Mito_percent > 10
# # table(discardCl)
# table((discardCl|discardDim|discardMito))
# seu_filt <- seu[, !(discardCl|discardDim|discardMito)]
# seu      <- seu_filt
# # seu_filt <- seu

# Save the processed Seurat object again after any optional manual filtering
saveRDS(seu, file = "results/01_Len_Integration.rds")

# Optional: Save session information for reproducibility
writeLines(capture.output(sessionInfo()), "session_info.txt")
