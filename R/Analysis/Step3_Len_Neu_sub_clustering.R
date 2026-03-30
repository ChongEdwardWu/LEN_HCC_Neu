### ==========================================================================
### Step 3: Further clustering and annotation on neutrophils
### ==========================================================================

### ---------------- Section 0: Preparation / Global options -----------------

# WARNING: this removes all existing objects in the R session
rm(list = ls())
gc()

# Optional: load a personal R profile if you use one locally
if (file.exists("path_to_.radian_profile")) {
  suppressMessages(source("path_to_.radian_profile"))
}

# ---------------------------------------------------------------------------
# Working directory and output folders
# ---------------------------------------------------------------------------
workdir <- "path_to_data"
setwd(workdir)

if (!file.exists(file.path(workdir, "figures"))) {
  dir.create(file.path(workdir, "figures"))
}
if (!file.exists(file.path(workdir, "results"))) {
  dir.create(file.path(workdir, "results"))
}

# ---------------------------------------------------------------------------
# Packages
# ---------------------------------------------------------------------------
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
  library(UCell)
  library(SingleCellExperiment)
  library(destiny)
  library(msigdbr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(dplyr)
  library(stringr)
})

# ! httpgd::hgd()

# ---------------------------------------------------------------------------
# Parallelization and global options
# ---------------------------------------------------------------------------
set.seed(123)

nworkers <- 20
# plan("sequential")
# plan("multicore", workers = nworkers)

# Species flag (used e.g. for blacklist; "mm" for mouse, "hs" for human)
species <- "mm"   # choices: "mm" or "hs"


## =======================================================================
## Section 1: Load objects and subset neutrophils from our dataset
## =======================================================================

# 1.1 Load the annotated Seurat object from Step 2
seu <- readRDS("results/02_Annotation_LSK.rds")
message("Loaded integrated object: cells = ", ncol(seu),
        "  genes = ", nrow(seu))

# 1.2 Load Ng reference neutrophil object (processed from GSE243466)
seu_Ng <- readRDS(
  "path_to_ng_reference_rds"
)

message("Loaded reference object: cells = ", ncol(seu_Ng),
        "  genes = ", nrow(seu_Ng))

# 1.3 Keep a copy of the full object for safety
seu_full <- seu

# 1.4 Subset our object to neutrophils only based on broad annotation
#     (CellType_l1 was defined in Step2)
seu_neu <- subset(seu_full, subset = CellType_l1 == "Neu")

# Sanity checks
table(seu_neu$CellType_l1)
table(seu_neu$CellType_l2)
seu_neu$CellType_l1 <- factor(seu_neu$CellType_l1)
seu_neu$CellType_l2 <- factor(seu_neu$CellType_l2)

# Work with RNA assay for mapping (log-normalized workflow)
DefaultAssay(seu_Ng)  <- "RNA"
DefaultAssay(seu_neu) <- "RNA"


## =======================================================================
## Section 2: Recompute PCA/UMAP on Ng reference (log-normalization)
##            following Ng et al. (for MapQuery)
## =======================================================================

# 2.1 Normalize Ng reference (LogNormalize)
seu_Ng <- NormalizeData(
  seu_Ng,
  normalization.method = "LogNormalize",
  scale.factor         = 10000,
  verbose              = FALSE
)

# 2.2 Identify highly variable genes (HVGs)
seu_Ng <- FindVariableFeatures(
  seu_Ng,
  selection.method = "vst",
  nfeatures        = 2000,
  verbose          = FALSE
)

# 2.3 Scale data on HVGs
seu_Ng <- ScaleData(
  seu_Ng,
  features = VariableFeatures(seu_Ng),
  verbose  = FALSE
)

# 2.4 PCA on HVGs
seu_Ng <- RunPCA(
  seu_Ng,
  features = VariableFeatures(seu_Ng),
  npcs     = 50,
  verbose  = FALSE
)

# 2.5 Neighbor graph and clustering using dims 1:35 (as in the paper)
dims_ref_use <- 1:35

seu_Ng <- FindNeighbors(
  seu_Ng,
  reduction = "pca",
  dims      = dims_ref_use
)

seu_Ng <- FindClusters(
  seu_Ng,
  resolution = 0.8
)

# 2.6 UMAP on PCA with dims 1:35 and n.neighbors = 30
#     IMPORTANT: return.model = TRUE is needed for MapQuery.
seu_Ng <- RunUMAP(
  seu_Ng,
  reduction    = "pca",
  dims         = dims_ref_use,
  n.neighbors  = 30,
  min.dist     = 0.3,
  return.model = TRUE,
  seed.use     = 123
)

# Use ManuscriptClusters from the original metadata as the reference stage
Idents(seu_Ng) <- factor(
  seu_Ng$ManuscriptClusters,
  levels = c("preNeu", "IMM 1", "IMM 2",
             "MAT 1", "MAT 2", "MAT 3", "MAT 4", "MAT 5",
             "T1", "T2", "T3")
)

# Quick check of UMAP colored by ManuscriptClusters
DimPlot(
  seu_Ng,
  reduction = "umap",
  group.by  = "ManuscriptClusters",
  label     = TRUE
) + coord_fixed(1)


## =======================================================================
## Section 3: Seurat transfer mapping (integration-style)
##            SelectIntegrationFeatures + FindTransferAnchors + MapQuery
## =======================================================================

# Goal:
# - Use Ng neutrophils as reference (seu_Ng)
# - Map our neutrophils (seu_neu) into the Ng PCA/UMAP space
# - Transfer Ng ManuscriptClusters as a predicted maturation stage

# 3.1 Log-normalize and find HVGs on our neutrophils
seu_neu <- NormalizeData(
  seu_neu,
  normalization.method = "LogNormalize",
  scale.factor         = 10000,
  verbose              = FALSE
)

seu_neu <- FindVariableFeatures(
  seu_neu,
  selection.method = "vst",
  nfeatures        = 2000,
  verbose          = FALSE
)

seu_neu <- ScaleData(
  seu_neu,
  features = VariableFeatures(seu_neu),
  verbose  = FALSE
)

# 3.2 PCA on HVGs (we will later use dims 1:35 for anchors)
seu_neu <- RunPCA(
  seu_neu,
  features = VariableFeatures(seu_neu),
  npcs     = 50,
  verbose  = FALSE
)

# 3.3 Select shared features for anchor finding
features_anchor <- SelectIntegrationFeatures(
  object.list = list(seu_Ng, seu_neu),
  nfeatures   = 2000
)

# 3.4 Define common PCA dimensions to use (consistent with Ng's dims 1:35)
dims_use <- 1:35

# 3.5 Find transfer anchors (LogNormalize-mode)
anchors <- FindTransferAnchors(
  reference           = seu_Ng,
  query               = seu_neu,
  normalization.method= "LogNormalize",
  reference.reduction = "pca",
  dims                = dims_use,
  features            = features_anchor
)

# 3.6 MapQuery: project our neutrophils onto Ng UMAP and transfer ManuscriptClusters
seu_neu_mapped <- MapQuery(
  anchorset           = anchors,
  reference           = seu_Ng,
  query               = seu_neu,
  refdata             = list(Ng_stage = "ManuscriptClusters"),
  reference.reduction = "pca",
  reduction.model     = "umap"   # use the UMAP model stored in seu_Ng
)

# Copy predicted Ng stage into a convenient column
seu_neu_mapped$Ng_stage <- seu_neu_mapped$predicted.Ng_stage

# Sanity check: distribution of predicted Ng stages
table(seu_neu_mapped$Ng_stage)

# Quick visualization: Ng reference vs our mapped neutrophils on the same UMAP
p_ref <- DimPlot(
  seu_Ng,
  reduction = "umap",
  group.by  = "ManuscriptClusters",
  label     = TRUE
) + coord_fixed(1) +
  ggtitle("Ng reference neutrophils (ManuscriptClusters)")

p_qry <- DimPlot(
  seu_neu_mapped,
  reduction = "ref.umap",
  group.by  = "Ng_stage",
  label     = TRUE
) + coord_fixed(1) +
  ggtitle("Our neutrophils mapped onto Ng UMAP (predicted Ng_stage)")

p_ref + p_qry


## =======================================================================
## Section 4: Diffusion map on our neutrophils (own diffusion space)
## =======================================================================

# NOTE:
# - Here we compute a diffusion map only on our neutrophils, without
#   forcing them into the Ng diffusion space.
# - DC coordinates can then be interpreted using Ng_stage and the
#   maturation score (computed in Section 5).

DefaultAssay(seu_neu_mapped) <- "RNA"

# 4.1 Use current HVGs for diffusion map
hvgs_neu <- VariableFeatures(seu_neu_mapped)
length(hvgs_neu)

# 4.2 Convert to SingleCellExperiment for destiny
sce_neu <- as.SingleCellExperiment(seu_neu_mapped, assay = "RNA")
sce_neu_hv <- sce_neu[rownames(sce_neu) %in% hvgs_neu, ]

# 4.3 Run diffusion map (destiny)
dm_neu <- DiffusionMap(
  sce_neu_hv,
  sigma  = "local",
  k      = 50,
  n_eigs = 20
)

# 4.4 Extract DC1 and DC2
dm_neu_evecs <- eigenvectors(dm_neu)[, 1:2]
colnames(dm_neu_evecs) <- c("DC1", "DC2")

# 4.5 Attach diffusion map as a DimReduc object
seu_neu_mapped[["diffmap"]] <- CreateDimReducObject(
  embeddings = dm_neu_evecs,
  key       = "DC_",
  assay     = DefaultAssay(seu_neu_mapped)
)


## =======================================================================
## Section 5: Compute neutrophil maturation score (Table S1 + UCell)
##            (same signature as Ng; applied to seu_neu_mapped)
## =======================================================================

# Path to Table S1 Excel file
table_s1_file <- "path_to_ng_table_s1_xlsx"

# 5.1 Load Table S1 and extract gene symbols
table_s1 <- read_excel(
  path  = table_s1_file,
  sheet = 1
)
colnames(table_s1) <- c("Number", "GeneID")

maturation_genes_raw <- unique(table_s1$GeneID)
maturation_genes_raw <- maturation_genes_raw[!is.na(maturation_genes_raw)]

# 5.2 Case-insensitive mapping of gene symbols to our object (seu_neu_mapped)
genes_obj <- rownames(seu_neu_mapped)
match_idx <- match(toupper(maturation_genes_raw), toupper(genes_obj))

maturation_genes <- genes_obj[match_idx[!is.na(match_idx)]]
missing_genes    <- maturation_genes_raw[is.na(match_idx)]

cat("Table S1 genes:", length(maturation_genes_raw), "\n")
cat("Mapped in our dataset:", length(maturation_genes), "\n")
cat("Missing (ignored):", paste(missing_genes, collapse = ", "), "\n\n")

maturation_geneSets <- list(NeutrophilMaturation = maturation_genes)

# 5.3 Compute UCell score directly on seu_neu_mapped
seu_neu_mapped <- UCell::AddModuleScore_UCell(
  seu_neu_mapped,
  features = maturation_geneSets,
  assay    = "RNA",
  ncores   = nworkers
)

# UCell adds a column 'NeutrophilMaturation_UCell' to meta.data
summary(seu_neu_mapped$NeutrophilMaturation_UCell)


## =======================================================================
## Section 6: Visualizations (Ng UMAP + diffusion map + maturation)
## =======================================================================

# 6.1 Our neutrophils on Ng UMAP, colored by predicted Ng_stage
p_umap_Ng_stage <- DimPlot(
  seu_neu_mapped,
  reduction = "ref.umap",
  group.by  = "Ng_stage",
  label     = TRUE
) + coord_fixed(1) +
  ggtitle("Our neutrophils on Ng UMAP (predicted Ng_stage)")

# 6.2 Same UMAP, split by treatment group
p_umap_split <- DimPlot(
  seu_neu_mapped,
  reduction = "ref.umap",
  group.by  = "Ng_stage",
  split.by  = "group",
  label     = FALSE
) + coord_fixed(1) +
  ggtitle("Ng_stage across treatment groups on Ng UMAP")

# 6.3 Diffusion map colored by Ng_stage
p_dm_stage <- DimPlot(
  seu_neu_mapped,
  reduction = "diffmap",
  group.by  = "Ng_stage",
  label     = TRUE
) + coord_fixed(1) +
  ggtitle("Diffusion map of our neutrophils (colored by Ng_stage)")

# 6.4 Diffusion map colored by treatment group
p_dm_group <- DimPlot(
  seu_neu_mapped,
  reduction = "diffmap",
  group.by  = "group",
  label     = FALSE
) + coord_fixed(1) +
  ggtitle("Diffusion map of our neutrophils (colored by treatment group)")

DimPlot(
  seu_neu_mapped,
  reduction = "diffmap",
  group.by  = "Ng_stage",
  split.by  = "group",
  label     = FALSE
) + coord_fixed(1) +
  ggtitle("Diffusion map of our neutrophils (colored by treatment group)")



# 6.5 Maturation score vs DC1 (trajectory-like view)
df_dc <- FetchData(
  seu_neu_mapped,
  vars = c("DC_1", "DC_2", "NeutrophilMaturation_UCell", "group", "Ng_stage")
)

p_dc_maturation <- ggplot(
  df_dc,
  aes(x = DC_1, y = NeutrophilMaturation_UCell, color = group)
) +
  geom_point(alpha = 0.5, size = 0.5) +
  theme_bw() +
  xlab("Diffusion component 1") +
  ylab("Neutrophil maturation score (UCell)") +
  ggtitle("Maturation score along DC1 in our neutrophils")

# Print / save plots as needed
p_umap_Ng_stage
p_umap_split
p_dm_stage
p_dm_group
p_dc_maturation


## =======================================================================
## Section 7: Build FigS2A gene sets (MOESM4) + UCell scores
##            (compute scores first; NO auto plotting)
## =======================================================================

# ---------------------------- User inputs --------------------------------
xlsx_path <- "path_to_figs2a_moesm4_xlsx"
fig_tag_pattern <- "(Sup\\s*Fig\\s*2A|Fig\\.?\\s*S2A|FigS2A|S2A)"  # robust match

assay_use      <- "RNA"     # use RNA assay for signature scoring
prefix_S2A     <- "S2A_"    # prefix for UCell meta columns
min_genes_S2A  <- 10        # drop signatures with too few genes after intersect
max_genes_S2A  <- 5000      # optional upper bound
ncores_ucell   <- 20

# Output directories for FigS2A
out_dir_S2A  <- file.path(workdir, "results",  "FigS2A")
fig_dir_S2A  <- file.path(workdir, "figures", "FigS2A")
if (!dir.exists(out_dir_S2A)) dir.create(out_dir_S2A, recursive = TRUE)
if (!dir.exists(fig_dir_S2A)) dir.create(fig_dir_S2A, recursive = TRUE)

# ----------------------- Helper: safe set names ---------------------------
safe_name <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

# ---------------- Helper: map gene symbols to object features (case-insensitive)
# Return the object's feature names (preserve exact casing in rownames(seu))
map_to_obj_features <- function(genes, obj_features) {
  genes <- unique(na.omit(as.character(genes)))
  if (!length(genes)) return(character(0))

  g_u <- toupper(genes)
  f_u <- toupper(obj_features)

  hit <- obj_features[match(g_u, f_u)]
  unique(na.omit(hit))
}

# ---------------- Helper: human -> mouse homolog mapping ------------------
# Prefer "homologene" if installed. If unavailable, keep original symbols.
map_hs_to_mm <- function(genes_hs) {
  genes_hs <- unique(na.omit(as.character(genes_hs)))

  if (requireNamespace("homologene", quietly = TRUE)) {
    hm <- homologene::homologene(genes_hs, inTax = 9606, outTax = 10090)
    out <- unique(na.omit(hm$mouseGene))
    if (!length(out)) warning("Homolog mapping returned 0 genes; keeping original symbols.")
    return(if (length(out)) out else genes_hs)
  }

  warning("Package 'homologene' not installed; returning human symbols unchanged.")
  genes_hs
}

# ---------------- Helper: resolve MSigDB gs_name candidates ---------------
make_msig_candidates <- function(gene_set_id) {
  x <- as.character(gene_set_id)
  x <- stringr::str_replace_all(x, "[\\.\\-\\s/]+", "_")
  x <- toupper(x)

  cand <- c(x)

  # Common pattern in MOESM4: "GSEA_XXX" should often map to MSigDB names without GSEA_
  if (stringr::str_detect(x, "^GSEA_")) {
    x2 <- stringr::str_replace(x, "^GSEA_", "")
    cand <- c(cand, x2)

    # Try GO BP prefix
    cand <- c(cand, paste0("GOBP_", x2))

    # Also try converting GSEA_ -> GOBP_ directly
    cand <- c(cand, stringr::str_replace(x, "^GSEA_", "GOBP_"))

    # Hallmark case: GSEA_HALLMARK_...
    cand <- c(cand, stringr::str_replace(x, "^GSEA_HALLMARK_", "HALLMARK_"))
  }

  # GO-like IDs may be "GO_..." -> msigdbr uses "GOBP_..."
  if (stringr::str_detect(x, "^GO_")) {
    cand <- c(cand, stringr::str_replace(x, "^GO_", "GOBP_"))
  }

  # Plain term -> also try "GOBP_"
  if (!stringr::str_detect(x, "^(HALLMARK|GOBP|REACTOME|KEGG|WP|PID|BIOCARTA)_")) {
    cand <- c(cand, paste0("GOBP_", x))
  }

  unique(cand)
}

# ---------------- Helper: fetch genes from a GO ID (fallback) -------------
fetch_go_genes_mouse <- function(go_id) {
  go_id <- as.character(go_id)

  tmp <- suppressMessages(suppressWarnings(
    AnnotationDbi::select(
      org.Mm.eg.db,
      keys    = go_id,
      keytype = "GOALL",
      columns = c("SYMBOL")
    )
  ))

  tmp %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::distinct(SYMBOL) %>%
    dplyr::pull(SYMBOL)
}


# -------------------------- 7.1 Read tables -------------------------------
# Read MOESM4 sheets using openxlsx2
gs_col <- openxlsx2::read_xlsx(xlsx_path, sheet = "Gene Sets collections")
gs_pub <- openxlsx2::read_xlsx(xlsx_path, sheet = "Gene Sets from publications")

# Robustly rename the "Shown in" column (openxlsx2 may trim trailing spaces)
shown_col <- grep("^Shown", colnames(gs_col), value = TRUE)
stopifnot(length(shown_col) == 1)

gs_col <- gs_col %>%
  dplyr::rename(
    shown_in = all_of(shown_col),
    species  = `Sp.`
  )

# Quick check (optional)
colnames(gs_col)

# Filter for FigS2A
col_S2A <- gs_col %>%
  dplyr::filter(!is.na(shown_in)) %>%
  dplyr::filter(stringr::str_detect(shown_in, stringr::regex(fig_tag_pattern, ignore_case = TRUE)))

pub_S2A <- gs_pub %>%
  dplyr::filter(!is.na(Notes)) %>%
  dplyr::filter(stringr::str_detect(Notes, stringr::regex(fig_tag_pattern, ignore_case = TRUE)))

message("Collections (FigS2A): ", nrow(col_S2A))
message("Publications rows (FigS2A): ", nrow(pub_S2A),
        " ; unique lists: ", dplyr::n_distinct(pub_S2A$list))


# -------------------------- 7.2 Build collection gene sets ----------------
# Collect msigdbr rows needed (mouse GO:BP; human GO:BP; human Hallmark)
m_mouse_go  <- msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
m_human_go  <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
m_human_h   <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")

ms_mouse_go <- m_mouse_go %>%
  dplyr::distinct(gs_name, gene_symbol)

ms_human <- dplyr::bind_rows(m_human_go, m_human_h) %>%
  dplyr::distinct(gs_name, gene_symbol)

gene_sets_col <- list()
meta_col <- tibble::tibble()

for (i in seq_len(nrow(col_S2A))) {
  bp    <- col_S2A$`Biological Process`[i]
  gs_id <- col_S2A$`GENE SET ID`[i]
  sp    <- tolower(col_S2A$species[i])
  urlid <- col_S2A$`URL/ ID`[i]

  # Fix typo "Glycolisis" -> "Glycolysis"
  bp_fixed <- ifelse(bp == "Glycolisis", "Glycolysis", bp)
  set_name <- safe_name(bp_fixed)

  candidates <- make_msig_candidates(gs_id)

  # Choose msigdbr table based on declared species
  # If declared "human" but no hit in human msigdbr, we will try mouse msigdbr as fallback.
  ms_df_primary   <- if (sp == "mouse") ms_mouse_go else ms_human
  ms_df_fallback  <- if (sp == "mouse") NULL       else ms_mouse_go

  # Find the first matching gs_name in primary
ms_names_u <- toupper(unique(ms_df_primary$gs_name))
hit_u <- intersect(toupper(candidates), ms_names_u)[1]

genes <- character(0)
source_used <- NA_character_

if (!is.na(hit_u)) {
  hit <- unique(ms_df_primary$gs_name)[match(hit_u, ms_names_u)]
  genes <- ms_df_primary %>%
    dplyr::filter(gs_name == hit) %>%
    dplyr::pull(gene_symbol) %>%
    unique() %>% na.omit()
  source_used <- paste0("msigdbr:", hit, " (", sp, ")")
} else {
  # If human msigdbr failed, try mouse msigdbr directly (better than 0 homologs)
  if (!is.null(ms_df_fallback)) {
    ms_names_u2 <- toupper(unique(ms_df_fallback$gs_name))
    hit_u2 <- intersect(toupper(candidates), ms_names_u2)[1]
    if (!is.na(hit_u2)) {
      hit2 <- unique(ms_df_fallback$gs_name)[match(hit_u2, ms_names_u2)]
      genes <- ms_df_fallback %>%
        dplyr::filter(gs_name == hit2) %>%
        dplyr::pull(gene_symbol) %>%
        unique() %>% na.omit()
      source_used <- paste0("msigdbr:", hit2, " (mouse fallback)")
    }
  }

  # Still no hit -> GO fallback if available
  if (!length(genes)) {
    if (is.character(urlid) && stringr::str_detect(urlid, "^GO:")) {
      genes <- fetch_go_genes_mouse(urlid)
      source_used <- paste0(urlid, " (org.Mm.eg.db fallback)")
    } else {
      warning("Could not resolve set from msigdbr and no GO ID fallback: ",
              bp_fixed, " | ", gs_id)
    }
  }
}


  # If declared as Human, map to mouse homologs
  if (length(genes) && sp == "human" && stringr::str_detect(source_used, "\\(human\\)")) {
  genes <- map_hs_to_mm(genes)
}


  gene_sets_col[[set_name]] <- unique(na.omit(genes))

  meta_col <- dplyr::bind_rows(meta_col, tibble::tibble(
    gene_set         = set_name,
    label_in_sheet   = bp_fixed,
    origin_sheet     = "Gene Sets collections",
    declared_species = sp,
    source           = source_used,
    n_genes          = length(gene_sets_col[[set_name]])
  ))
}


# -------------------------- 7.3 Build publication gene sets ---------------
# Make sure gs_pub has a "species" column too (optional but consistent)
# If your publications sheet already has "Sp.", rename it here; otherwise ignore.
if ("Sp." %in% colnames(pub_S2A) && !"species" %in% colnames(pub_S2A)) {
  pub_S2A <- pub_S2A %>% dplyr::rename(species = `Sp.`)
}

pub_tbl <- pub_S2A %>%
  dplyr::group_by(list) %>%
  dplyr::summarise(
    genes            = list(unique(na.omit(gene))),
    declared_species = tolower(dplyr::first(species)),
    source           = paste(unique(source), collapse = "; "),
    .groups          = "drop"
  ) %>%
  dplyr::mutate(gene_set = safe_name(list))

gene_sets_pub <- stats::setNames(pub_tbl$genes, pub_tbl$gene_set)

# Map any human publication list to mouse homologs
is_human_pub <- pub_tbl$declared_species == "human"
if (any(is_human_pub)) {
  for (nm in pub_tbl$gene_set[is_human_pub]) {
    gene_sets_pub[[nm]] <- map_hs_to_mm(gene_sets_pub[[nm]])
  }
}

meta_pub <- pub_tbl %>%
  dplyr::transmute(
    gene_set,
    label_in_sheet   = list,
    origin_sheet     = "Gene Sets from publications",
    declared_species,
    source,
    n_genes          = lengths(gene_sets_pub[gene_set])
  )


# -------------------------- 7.4 Merge + export (raw) ----------------------
gene_sets_S2A_raw <- c(gene_sets_col, gene_sets_pub)

# Avoid rare name collisions
names(gene_sets_S2A_raw) <- make.unique(names(gene_sets_S2A_raw))

meta_S2A <- dplyr::bind_rows(meta_col, meta_pub) %>%
  dplyr::arrange(dplyr::desc(n_genes))

gene_sets_long <- tibble::enframe(gene_sets_S2A_raw, name = "gene_set", value = "gene") %>%
  tidyr::unnest(gene)

saveRDS(gene_sets_S2A_raw, file.path(out_dir_S2A, "FigS2A_gene_sets_list_raw.rds"))
utils::write.csv(meta_S2A,  file.path(out_dir_S2A, "FigS2A_gene_sets_summary.csv"),
                 row.names = FALSE)
utils::write.csv(gene_sets_long, file.path(out_dir_S2A, "FigS2A_gene_sets_long.csv"),
                 row.names = FALSE)

message("FigS2A gene sets (raw) built: ", length(gene_sets_S2A_raw))


# -------------------------- 7.5 UCell scoring (NO plotting) ---------------
obj_name <- "seu_neu_mapped"
seu_tmp <- seu_neu_mapped

assay_sc <- assay_use
if (!assay_sc %in% names(seu_tmp@assays)) assay_sc <- "RNA"
Seurat::DefaultAssay(seu_tmp) <- assay_sc

# remove old S2A columns if re-run
old_cols <- grep(paste0("^", prefix_S2A, ".*_UCell$"),
                 colnames(seu_tmp@meta.data), value = TRUE)
if (length(old_cols) > 0) seu_tmp@meta.data[, old_cols] <- NULL

obj_features <- rownames(seu_tmp)

# case-insensitive intersect to object features
gs_in_obj <- lapply(gene_sets_S2A_raw, map_to_obj_features, obj_features = obj_features)
n_in_obj  <- lengths(gs_in_obj)

qc_tbl <- tibble::tibble(
  gene_set          = names(gene_sets_S2A_raw),
  n_genes_raw       = lengths(gene_sets_S2A_raw),
  n_genes_in_object = n_in_obj
) %>%
  dplyr::arrange(dplyr::desc(n_genes_in_object), dplyr::desc(n_genes_raw))

utils::write.csv(qc_tbl,
                 file.path(out_dir_S2A, "FigS2A_gene_sets_QC_in_object_seu_neu_mapped.csv"),
                 row.names = FALSE)

keep <- (n_in_obj >= min_genes_S2A) & (n_in_obj <= max_genes_S2A)
gs_use <- gs_in_obj[keep]

saveRDS(gs_use, file.path(out_dir_S2A, "FigS2A_gene_sets_list_USED_in_seu_neu_mapped.rds"))

message("Gene sets retained for UCell = ", length(gs_use))

gs_ucell <- gs_use
names(gs_ucell) <- paste0(prefix_S2A, names(gs_ucell))

seu_tmp <- UCell::AddModuleScore_UCell(
  seu_tmp,
  features = gs_ucell,
  assay    = assay_sc,
  ncores   = ncores_ucell
)

# assign back + save
seu_neu_mapped <- seu_tmp
rm(seu_tmp)

out_rds <- file.path(out_dir_S2A, "seu_neu_mapped_with_FigS2A_UCell.rds")
saveRDS(seu_neu_mapped, out_rds)
message("Saved scored object -> ", out_rds)


## =======================================================================
## Section 8: Plot FigS2A-like panel on seu_neu_mapped
## =======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(scales)
})

reduction_use <- "ref.umap"   # or "umap"
pt_size_use   <- 0.15
qmin_use      <- "q05"
qmax_use      <- "q95"
ncol_use      <- 5

# Your preferred color scale
s2a_cols <- c('#1f618d', '#2E86C1', "white", '#ec7063', "#AD272B")

# Helper: one FeaturePlot with styling (NO duplicate scale warning)
plot_one <- function(seu, title, feature, reduction = "ref.umap",
                     pt.size = 0.15, qmin = "q05", qmax = "q95") {

  md <- colnames(seu@meta.data)
  if (!feature %in% md) {
    message("Skip (missing): ", title, " | ", feature)
    return(NULL)
  }

  p <- Seurat::FeaturePlot(
    seu,
    reduction  = reduction,
    features   = feature,
    pt.size    = pt.size,
    order      = TRUE,
    min.cutoff = qmin,
    max.cutoff = qmax
  )

  # Remove Seurat's default scale and add yours (avoid "Scale already present" messages)
  p <- p +
    ggplot2::scale_color_gradientn(colours = s2a_cols) +
    coord_fixed(1) +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9),
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank()
    ) +
    NoLegend()

  p
}

# --------------------- Panel order (use your exact column names) ----------
panel_tbl <- tibble::tribble(
  ~title,                ~feature,
  "mRNA counts",         "nCount_RNA",
  "ROS production",      "S2A_ROS_production_UCell",
  "Glycolysis",          "S2A_Glycolysis_UCell",
  "Response to lipids",  "S2A_Response_to_lipids_UCell",
  "IFN-I response",      "S2A_TypeI_IFN_UCell",
  "Extravasation",       "S2A_Extravasation_UCell",
  "Chronic Inflamm.",    "S2A_Chronic_Inflammation_UCell",
  "NETosis",             "S2A_NETosis_UCell",
  "DNA damage",          "S2A_DNA_damage_UCell",
  "Cell death",          "S2A_Progr_cell_death_UCell",

  "Phagocytosis",        "S2A_Phagocytosis_UCell",
  "Cell-cell adhesion",  "S2A_Adhesion_UCell",
  "Autoimmunity",        "S2A_Autoimmunity_UCell",
  "PMN-MDSC",            "S2A_PMN_MDSC_UCell",
  "Activated\nPMN-MDSC", "S2A_Activated_PMN_MDSC_UCell",
  "TAN Lung",            "S2A_TANS_LUNG_UCell",
  "PDAC (TAN1)",         "S2A_PDAC_TAN_1_UCell",
  "APCTAN",              "S2A_APC_TANs_UCell",
  "PDL1",                "S2A_PDL1_NEUTROPHILS_UCell"
)

# Build plots
plots <- purrr::pmap(
  panel_tbl,
  ~ plot_one(seu_neu_mapped, title = ..1, feature = ..2,
             reduction = reduction_use, pt.size = pt_size_use,
             qmin = qmin_use, qmax = qmax_use)
)

plots <- Filter(Negate(is.null), plots)
p_S2A_like <- patchwork::wrap_plots(plots, ncol = ncol_use)
p_S2A_like

# Save
ggsave(
  filename = file.path(fig_dir_S2A, paste0("FigS2A_like_projection_", reduction_use, ".png")),
  plot     = p_S2A_like,
  width    = 16, height = 7, dpi = 300
)
ggsave(
  filename = file.path(fig_dir_S2A, paste0("FigS2A_like_projection_", reduction_use, ".pdf")),
  plot     = p_S2A_like,
  width    = 16, height = 7
)

## =======================================================================
## Section 9: Save objects and session information
## =======================================================================

# Save the mapped neutrophil object
saveRDS(seu_neu_mapped, file = "results/03_Len_Neu_NgMapping.rds")
# seu_neu_mapped <- readRDS("results/03_Len_Neu_NgMapping.rds")

# Optionally save the updated Ng reference (with new PCA/UMAP)
saveRDS(seu_Ng, file = "results/03_Ng_ref_updated.rds")
# seu_Ng <- readRDS("results/03_Ng_ref_updated.rds")

# Save session information for reproducibility
writeLines(capture.output(sessionInfo()), "session_info.txt")
