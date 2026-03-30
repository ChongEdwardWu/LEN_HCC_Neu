### Section 0: Preparation ----------------------------------------------------
rm(list = ls()); gc()

# Optional: load a personal R profile if you use one locally
if (file.exists("path_to_.radian_profile")) {
  suppressMessages(source("path_to_.radian_profile"))
}
# .libPaths()

# Working directory
workdir <- "path_to_data"
setwd(workdir)
if (!file.exists(file.path(workdir, "figures"))) dir.create(file.path(workdir, "figures"))
if (!file.exists(file.path(workdir, "results"))) dir.create(file.path(workdir, "results"))

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
  library(ggrepel)
  # Optional blocks:
  # library(GSEABase); library(AUCell)
  # library(clusterProfiler); library(org.Mm.eg.db)
})

set.seed(123)
nworkers <- 8

## Species flag
species <- "mm"
# species <- "hs"


# Create directories for figures and results
if (file.exists(file.path(workdir, "figures/final"))) {
} else {
  dir.create(file.path(workdir, "figures/final"))
}

# Load Seurat objects used by the figure panels
seu_full <- readRDS("results/02_Annotation_LSK.rds")
seu_neu <- readRDS("results/03_Len_Neu_NgMapping.rds")

# ! httpgd::hgd()


### Section 1: Figure 1. All cells dimension reduction ------------------------
seu <- seu_full[,seu_full$group %in% c("CTR", "LEN")]

table(seu$CellType_l1)

Idents(seu) <- seu$CellType_l1
# define colors
colorpanel <- c(
    "T" = "#575bb7",
    "B" = "#8C7A61",
    "NK" = "#4FA5F6",
    "DC" = "#C4DA5D", # lime-chartreuse
    "pDC" = "#54a546",
    "Mono" = "#66cccc",
    "Mac" = "#fa9645",
    "Neu" = "#e13d2d",
    "Baso" = "#E053B2" # magenta
)


# Plot 1: Umap and clusters 
DimPlot(seu,
        reduction = "umap",
        # split.by = "group",
        group.by = "CellType_l1", 
        label = F
) +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = colorpanel)

ggsave(file = paste0("figures/final/fig1_all_UMAP_L1.png"), width = 150, height = 150, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/final/fig1_all_UMAP_L1.pdf"), width = 150, height = 150, units = "mm", device = "pdf", bg = 'transparent')


### Section 2: Figure 2. Neutrophil dimension reduction -----------------------
seu <- seu_neu[,seu_neu$group %in% c("CTR", "LEN")]

table(seu$Ng_stage)

# Recode Ng_stage -> CellType_l2
seu$CellType_l2 <- dplyr::recode(
  as.character(seu$Ng_stage),
  `preNeu` = "preNeu",
  `IMM 1`  = "IMM",
  `IMM 2`  = "IMM",
  `MAT 1`  = "MAT",
  `MAT 2`  = "MAT",
  `MAT 3`  = "MAT",
  `MAT 4` = "MAT",
  `MAT 5` = "MAT",
  `T1`     = "IMM_TAN",
  `T2`     = "MAT_TAN",
  `T3`     = "IMM_TAN",
  .default = NA_character_
)

# Set factor levels order
seu$CellType_l2 <- factor(seu$CellType_l2, levels = c("preNeu", "IMM", "MAT", "IMM_TAN","MAT_TAN"))

# Quick check
table(seu$Ng_stage, seu$CellType_l2, useNA = "ifany")
table(seu$CellType_l2, useNA = "ifany")

Idents(seu) <- seu$CellType_l2

# define colors
colorpanel <- c(
    "preNeu" = "#5658c7",
    "IMM" = "#4fcaf6",
    "MAT" = "#2ccc8f",
    "IMM_TAN" = "#F09000",
    "MAT_TAN" = "#af3c13"
)


# Plot 2: Umap and clusters 
DimPlot(seu,
        reduction = "ref.umap",
        split.by = "group",
        group.by = "CellType_l2", 
        label = F
) +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = colorpanel)

ggsave(file = paste0("figures/final/fig2_neu_refUMAP_L2.png"), width = 250, height = 150, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/final/fig2_neu_refUMAP_L2.pdf"), width = 250, height = 150, units = "mm", device = "pdf", bg = 'transparent')


### Section 3: Figure 3. Cluster abundance -----------------------------------

# View and plot cluster distribution
source_cluster <- tibble(
  seu$CellType_l2,
  seu$group
) %>%
  set_names("Cluster", "Group") %>%
  group_by(Cluster, Group) %>%
  summarise(no.cell = n()) %>%
  group_by(Group) %>%
  mutate(
    total.no = sum(no.cell),
    perc = 100 * no.cell / total.no
  )
# View(source_cluster)

source_cluster_print <- source_cluster %>%
  ungroup() %>%
  mutate(perc = round(perc, 2)) %>%
  arrange(Group, desc(perc), Cluster)

cat("\nCellType_l2 proportions within each group (%):\n")
print(source_cluster_print, n = Inf)


library(ggplot2)
ggplot(source_cluster, aes(x = Group, y = perc, fill = Cluster)) +
  geom_col(colour = "black") +
  scale_fill_manual(values=colorpanel) +
  coord_fixed(ratio = 1 / 10) +
  theme_bw() +
  xlab("group") +
  ylab("%")


ggsave(file = paste0("figures/final/fig3_Neu_L2_cluster_abundance.png"), width = 100, height = 100, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/final/fig3_Neu_L2_cluster_abundance.pdf"), width = 100, height = 100, units = "mm", device = "pdf", bg = 'transparent')


### Section 4: Figure 4. Neutrophil maturation score ridge plot ----------------

p_ridge_mat <- Seurat::RidgePlot(
  seu,
  features = "NeutrophilMaturation_UCell",
  group.by = "CellType_l2",
  ncol = 1,
  same.y.lims = TRUE
) +
  ggtitle("Maturation signature (NeutrophilMaturation_UCell)") +
  xlab("Maturation score") +
  ylab(NULL) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.y = element_blank()
  ) +
  scale_fill_manual(values=colorpanel)

p_ridge_mat


ggsave(file = paste0("figures/final/fig4_Neu_L2_maturation_score.png"), width = 150, height = 120, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/final/fig4_Neu_L2_maturation_score.pdf"), width = 150, height = 120, units = "mm", device = "pdf", bg = 'transparent')


### Section 5: Figure 5. DEG comparison --------------------------------------
bulkDEG <- openxlsx2::read_xlsx("path_to_bulk_deg_xlsx", sheet = "DEG")
scDEG <- openxlsx2::read_xlsx(
  file.path(workdir, "results", "Step4_GroupDE_GSEA", "Step4_LEN_vs_CTR_AllNeu.xlsx"),
  sheet = "DEG_RNA"
)

# --------------------------- Inputs ---------------------------------------
padj_th_bulk <- 0.05
padj_th_sc   <- 0.05

# If you want comparable effect-size thresholds, set these:
lfc_th_bulk  <- 0.5   # or 0.25 if you want stricter
lfc_th_sc    <- 0.25   # sc already had 0.25 in FindMarkers output, but keep as 0 here

# Your dataframes: bulkDEG, scDEG
# bulkDEG columns assumed: GeneSymbol, log2FoldChange, padj
# scDEG columns assumed: feature, avg_log2FC, p_val_adj

# --------------------------- Tidy + merge ---------------------------------
bulk_tbl <- bulkDEG %>%
  dplyr::transmute(
    gene_bulk      = GeneSymbol,
    gene_key       = toupper(GeneSymbol),
    bulk_log2FC    = log2FoldChange,
    bulk_padj      = padj
  ) %>%
  dplyr::filter(!is.na(gene_key), gene_key != "") %>%
  dplyr::distinct(gene_key, .keep_all = TRUE)

sc_tbl <- scDEG %>%
  dplyr::transmute(
    gene_sc        = feature,
    gene_key       = toupper(feature),
    sc_log2FC      = avg_log2FC,
    sc_padj        = p_val_adj
  ) %>%
  dplyr::filter(!is.na(gene_key), gene_key != "") %>%
  dplyr::distinct(gene_key, .keep_all = TRUE)

merged <- dplyr::full_join(bulk_tbl, sc_tbl, by = "gene_key") %>%
  dplyr::mutate(
    gene = dplyr::coalesce(gene_bulk, gene_sc),

    # significance flags
    bulk_sig = !is.na(bulk_padj) & bulk_padj < padj_th_bulk & !is.na(bulk_log2FC) & abs(bulk_log2FC) > lfc_th_bulk,
    sc_sig   = !is.na(sc_padj)   & sc_padj   < padj_th_sc   & !is.na(sc_log2FC)   & abs(sc_log2FC)   > lfc_th_sc,

    # direction
    bulk_dir = dplyr::case_when(
      !is.na(bulk_log2FC) & bulk_log2FC > 0 ~ "UP",
      !is.na(bulk_log2FC) & bulk_log2FC < 0 ~ "DN",
      TRUE ~ NA_character_
    ),
    sc_dir = dplyr::case_when(
      !is.na(sc_log2FC) & sc_log2FC > 0 ~ "UP",
      !is.na(sc_log2FC) & sc_log2FC < 0 ~ "DN",
      TRUE ~ NA_character_
    ),

    # category
    category = dplyr::case_when(
      bulk_sig & sc_sig & bulk_dir == "UP" & sc_dir == "UP" ~ "Concordant_UP",
      bulk_sig & sc_sig & bulk_dir == "DN" & sc_dir == "DN" ~ "Concordant_DN",
      bulk_sig & sc_sig & bulk_dir != sc_dir               ~ "Discordant",
      bulk_sig & !sc_sig                                   ~ "Bulk_only",
      !bulk_sig & sc_sig                                   ~ "sc_only",
      TRUE                                                 ~ "Not_sig"
    )
  )

# Quick summary
merged %>%
  dplyr::count(category) %>%
  dplyr::arrange(dplyr::desc(n))

concord_up <- merged %>%
  dplyr::filter(category == "Concordant_UP") %>%
  dplyr::arrange(bulk_padj, sc_padj)

concord_dn <- merged %>%
  dplyr::filter(category == "Concordant_DN") %>%
  dplyr::arrange(bulk_padj, sc_padj)

discordant <- merged %>%
  dplyr::filter(category == "Discordant") %>%
  dplyr::arrange(bulk_padj, sc_padj)

bulk_only <- merged %>%
  dplyr::filter(category == "Bulk_only") %>%
  dplyr::arrange(bulk_padj)

sc_only <- merged %>%
  dplyr::filter(category == "sc_only") %>%
  dplyr::arrange(sc_padj)

# Write to Excel
out_xlsx <- file.path(workdir, "results", "step5_Bulk_scRNA_merged_DEG_LEN_vs_CON.xlsx")
dir.create(dirname(out_xlsx), recursive = TRUE, showWarnings = FALSE)

openxlsx2::write_xlsx(
  x = list(
    Summary = merged %>% dplyr::count(category),
    Concordant_UP = concord_up,
    Concordant_DN = concord_dn,
    Discordant    = discordant,
    Bulk_only     = bulk_only,
    sc_only       = sc_only,
    All_merged    = merged
  ),
  file = out_xlsx
)
out_xlsx


# only genes present in both
df_scatter <- merged %>%
  dplyr::filter(!is.na(bulk_log2FC), !is.na(sc_log2FC)) %>%
  # Show only genes with sufficiently large bulk changes.
  dplyr::filter(abs(bulk_log2FC) >= lfc_th_bulk) %>%
  dplyr::mutate(
    plot_cat = dplyr::case_when(
      category == "Concordant_UP" ~ "Concordant_UP",
      category == "Concordant_DN" ~ "Concordant_DN",
      # Other points in the upper-right quadrant after excluding concordant genes.
      bulk_log2FC > 0 & sc_log2FC > 0 ~ "RU_other",
      # Other points in the lower-left quadrant.
      bulk_log2FC < 0 & sc_log2FC < 0 ~ "LD_other",
      TRUE ~ "Other"
    ),
    plot_cat = factor(
      plot_cat,
      levels = c("Concordant_UP", "Concordant_DN", "RU_other", "LD_other", "Other")
    )
  )

message("Points shown after bulk threshold filter (|bulk| >= ", lfc_th_bulk, "): ", nrow(df_scatter))

# Correlation is calculated on the displayed points only.
cor_spearman <- suppressWarnings(stats::cor(df_scatter$bulk_log2FC, df_scatter$sc_log2FC, method = "spearman"))
cor_pearson  <- suppressWarnings(stats::cor(df_scatter$bulk_log2FC, df_scatter$sc_log2FC, method = "pearson"))
message("Shown points: Spearman r = ", round(cor_spearman, 3), " | Pearson r = ", round(cor_pearson, 3))

# Label only Concordant_UP and Concordant_DN points.
df_label <- df_scatter %>%
  dplyr::filter(plot_cat %in% c("Concordant_UP", "Concordant_DN")) %>%
  dplyr::mutate(label = gene)

# If too many concordant genes are labeled, keep only the top N by absolute effect size.
# top_n_each <- 30
# df_label <- df_label %>%
#   dplyr::group_by(plot_cat) %>%
#   dplyr::slice_max(order_by = abs(bulk_log2FC) + abs(sc_log2FC),
#                    n = top_n_each, with_ties = FALSE) %>%
#   dplyr::ungroup()

# Color palette for the scatter plot.
cols_scatter <- c(
  Concordant_UP = "#AD272B",  # dark red
  Concordant_DN = "#1f618d",  # dark blue
  RU_other      = "#F4A6B6",  # pink
  LD_other      = "#A6D4FA",  # light blue
  Other         = "#BDBDBD"   # grey
)

# Plot the merged bulk/scRNA DEG scatter.
p_scatter <- ggplot2::ggplot(df_scatter, ggplot2::aes(x = bulk_log2FC, y = sc_log2FC)) +
  # zero axes
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  # bulk effect-size threshold lines
  ggplot2::geom_vline(xintercept = c(-lfc_th_bulk, lfc_th_bulk),
                      linetype = "dotted", linewidth = 0.55) +
  ggplot2::geom_point(ggplot2::aes(color = plot_cat), size = 0.9, alpha = 0.75) +
  ggrepel::geom_text_repel(
    data = df_label,
    ggplot2::aes(label = label, color = plot_cat),
    size = 3,
    max.overlaps = Inf,      # set to 50 or 100 if labels become too dense
    box.padding = 0.25,
    point.padding = 0.15,
    segment.size = 0.2,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = cols_scatter, drop = FALSE, name = NULL) +
  ggplot2::labs(
    title = paste0(
      "LEN vs CON: bulk vs scRNA log2FC (|bulk| ≥ ", lfc_th_bulk,
      "; Spearman=", round(cor_spearman, 3), ")"
    ),
    x = "Bulk log2FC (LEN - CON)",
    y = "scRNA avg_log2FC (LEN - CTR)"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "right")

p_scatter

ggsave(file = paste0("figures/final/fig5_Neu_DEG.png"), width = 200, height = 150, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/final/fig5_Neu_DEG.pdf"), width = 200, height = 150, units = "mm", device = "pdf", bg = 'transparent')


### Section 6: Figure 6. Gene signature score plot ---------------------------
library(clusterProfiler)
library(org.Mm.eg.db)    # for mouse
library(stringr)
library(fgsea)
library(msigdbr)
library(UCell)



Neupathways <- readRDS(file.path(workdir, "results", "FigS2A", "FigS2A_gene_sets_list_raw.rds"))
pathways <- Neupathways[c("PMN_MDSC","Extravasation")]
pathwaysH <- msigdbr("mouse", category = "H")
pathwaysH <- split(pathwaysH$gene_symbol, pathwaysH$gs_name)
pathways <- c(pathways, pathwaysH["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
pathwaysGOBP <- msigdbr("mouse", category = "C5", subcategory = "GO:BP")
pathwaysGOBP <- split(pathwaysGOBP$gene_symbol, pathwaysGOBP$gs_name)
pathways <- c(pathways, pathwaysGOBP["GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"])

# UCell_cols <- grep("_UCell$", colnames(seu@meta.data), value=T)
# seu@meta.data[,UCell_cols] <- NULL


seu <- UCell::AddModuleScore_UCell(
  seu,
  features = pathways,
  assay    = "RNA",
  ncores   = 16
)



# ---------------------------- Inputs -------------------------------------
UCell_cols <- grep("_UCell$", colnames(seu@meta.data), value=T)

# Group order (edit if you later include all 6 groups)
group_levels <- c("CTR", "LEN")
group_levels <- intersect(group_levels, unique(as.character(seu$group)))
seu$group <- factor(as.character(seu$group), levels = group_levels)

# Pretty labels (optional)
sig_pretty <- c(
  PMN_MDSC_UCell = "PMN-MDSC",
  Extravasation_UCell   = "Extravasation",
  HALLMARK_INTERFERON_GAMMA_RESPONSE_UCell = "INTERFERON_GAMMA_RESPONSE",
  GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_UCell = "ANTIGEN_PROCESSING_AND_PRESENTATION"
)


df_long <- Seurat::FetchData(seu, vars = c("group", UCell_cols)) %>%
  tibble::as_tibble() %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(UCell_cols),
    names_to = "Signature",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    group = factor(as.character(group), levels = group_levels),
    Signature = factor(Signature, levels = UCell_cols),
    Signature_pretty = factor(
      dplyr::recode(as.character(Signature), !!!sig_pretty),
      levels = unname(sig_pretty[UCell_cols])
    )
  )

# Optional: set clearer group colors (edit if you want)
grp_fill <- c(CTR = "#A6D4FA", LEN = "#F4A6B6")

p_vln <- ggplot2::ggplot(df_long, ggplot2::aes(x = group, y = Score, fill = group)) +
  ggplot2::geom_violin(trim = FALSE, scale = "width", color = "black") +
  ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, color = "black") +
  ggplot2::stat_summary(fun = median, geom = "point", size = 1.2, color = "black") +
  ggplot2::facet_wrap(~ Signature_pretty, scales = "free_y", ncol = 2) +
  ggplot2::scale_fill_manual(values = grp_fill) +
  ggplot2::labs(x = NULL, y = "UCell score") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(face = "bold", size = 11),
    axis.text.x = ggplot2::element_text(hjust = 0.5)
  ) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.08)))

p_vln

ggsave(file = paste0("figures/final/fig6_Neu_UCell_VlnPlot.png"), p_vln, width = 250, height = 150, units = "mm", dpi = 300, device = "png")
ggsave(file = paste0("figures/final/fig6_Neu_UCell_VlnPlot.pdf"), p_vln, width = 250, height = 150, units = "mm", device = "pdf", bg = 'transparent')


## =======================================================================
### Section 7: Figures 7-8. Venn overlap and DEG export -----------------------
## Input: Step4_GroupDE_GSEA Excel files (DEG_RNA sheets)
## Output: Venn PNG/PDF + overlap gene lists + DEG detail tables


## Venn and heatmap setup -----------------------


suppressPackageStartupMessages({
  library(stringr)
  library(tibble)
  library(purrr)
  library(openxlsx2)
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(pheatmap)
  library(grid)
  library(VennDiagram)
})

## User parameters --------------------------------------------------------
in_dir <- file.path(workdir, "results", "Step4_GroupDE_GSEA")

xlsx_len <- file.path(in_dir, "Step4_LEN_vs_CTR_AllNeu.xlsx")
xlsx_amg <- file.path(in_dir, "Step4_LENnAMG_vs_LEN_AllNeu.xlsx")
xlsx_ko  <- file.path(in_dir, "Step4_LENnKO_vs_LEN_AllNeu.xlsx")

# Save final figures to figures/final (as requested)
fig_dir <- file.path(workdir, "figures", "final")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Results tables are saved under results/
out_dir <- file.path(workdir, "results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)



## DEG thresholds (adjust if you want)
padj_cut  <- 0.05
logfc_cut <- 0.25

## Heatmap settings
heat_groups <- c("LEN", "LENnKO", "LENnAMG")   # row order like paper: Len, Len&CKO, Len&AMG
heat_group_labels <- c("Len", "Len&CKO", "Len&AMG")

## Seurat object for heatmap expression
## Use the same RDS you used in Step4 (should contain all 6 groups).
seu_rds <- file.path(workdir, "results", "03_Len_Neu_NgMapping.rds")

## Helper functions -------------------------------------------------------

# 1.1 Read a sheet from Step4 excel (openxlsx2)
read_step4_sheet <- function(xlsx, sheet) {
  df <- openxlsx2::read_xlsx(xlsx, sheet = sheet)
  df <- tibble::as_tibble(df)
  df
}

# 1.2 Robustly detect logFC column name from Seurat FindMarkers outputs
detect_fc_col <- function(df) {
  cand <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  hit <- intersect(cand, colnames(df))
  if (length(hit) == 0) {
    stop("Cannot find a logFC column in DEG table. Columns are: ",
         paste(colnames(df), collapse = ", "))
  }
  hit[1]
}

# 1.3 Get significant up/down gene sets
get_deg_genes <- function(deg_df,
                          direction = c("up", "down"),
                          padj_cut = 0.05,
                          logfc_cut = 0.25,
                          feature_col = "feature") {
  direction <- match.arg(direction)

  if (!feature_col %in% colnames(deg_df)) {
    stop("DEG table does not contain column '", feature_col, "'. Columns: ",
         paste(colnames(deg_df), collapse = ", "))
  }
  fc_col <- detect_fc_col(deg_df)
  if (!"p_val_adj" %in% colnames(deg_df)) {
    stop("DEG table does not contain column 'p_val_adj'. Columns: ",
         paste(colnames(deg_df), collapse = ", "))
  }

  deg_df <- deg_df %>%
    dplyr::filter(!is.na(.data[[feature_col]])) %>%
    dplyr::filter(!is.na(p_val_adj)) %>%
    dplyr::filter(p_val_adj <= padj_cut) %>%
    dplyr::filter(!is.na(.data[[fc_col]]))

  if (direction == "up") {
    out <- deg_df %>%
      dplyr::filter(.data[[fc_col]] >= logfc_cut) %>%
      dplyr::pull(.data[[feature_col]]) %>%
      unique()
  } else {
    out <- deg_df %>%
      dplyr::filter(.data[[fc_col]] <= -logfc_cut) %>%
      dplyr::pull(.data[[feature_col]]) %>%
      unique()
  }
  out
}

# 1.4 Keep full DEG rows for the gene sets used in Venn diagrams
get_deg_table <- function(deg_df,
                          direction = c("up", "down"),
                          padj_cut = 0.05,
                          logfc_cut = 0.25,
                          feature_col = "feature",
                          comparison = NA_character_,
                          venn_panel = NA_character_,
                          set_label = NA_character_) {
  direction <- match.arg(direction)

  if (!feature_col %in% colnames(deg_df)) {
    stop("DEG table does not contain column '", feature_col, "'. Columns: ",
         paste(colnames(deg_df), collapse = ", "))
  }
  fc_col <- detect_fc_col(deg_df)
  if (!"p_val_adj" %in% colnames(deg_df)) {
    stop("DEG table does not contain column 'p_val_adj'. Columns: ",
         paste(colnames(deg_df), collapse = ", "))
  }

  out <- deg_df %>%
    dplyr::filter(!is.na(.data[[feature_col]])) %>%
    dplyr::filter(!is.na(p_val_adj)) %>%
    dplyr::filter(p_val_adj <= padj_cut) %>%
    dplyr::filter(!is.na(.data[[fc_col]]))

  if (direction == "up") {
    out <- out %>%
      dplyr::filter(.data[[fc_col]] >= logfc_cut)
  } else {
    out <- out %>%
      dplyr::filter(.data[[fc_col]] <= -logfc_cut)
  }

  out %>%
    dplyr::distinct(.data[[feature_col]], .keep_all = TRUE) %>%
    dplyr::mutate(
      comparison = comparison,
      direction = direction,
      venn_panel = venn_panel,
      set_label = set_label,
      padj_cutoff = padj_cut,
      logfc_cutoff = logfc_cut
    ) %>%
    dplyr::relocate(comparison, direction, venn_panel, set_label,
                    padj_cutoff, logfc_cutoff, dplyr::all_of(feature_col))
}

# 1.5 Build a Venn membership table so each gene can be traced back to a region
build_venn_membership <- function(sets_named_list, membership_labels) {
  if (length(sets_named_list) != length(membership_labels)) {
    stop("membership_labels must have the same length as sets_named_list.")
  }

  all_genes <- sort(unique(unlist(sets_named_list, use.names = FALSE)))
  if (length(all_genes) == 0) {
    return(tibble::tibble(gene = character(0)))
  }

  out <- tibble::tibble(gene = all_genes)
  for (i in seq_along(sets_named_list)) {
    out[[paste0("in_", membership_labels[[i]])]] <- all_genes %in% sets_named_list[[i]]
  }

  membership_cols <- paste0("in_", membership_labels)
  out$n_sets <- rowSums(as.matrix(out[, membership_cols, drop = FALSE]))
  out$region <- vapply(seq_len(nrow(out)), function(i) {
    hits <- membership_labels[as.logical(out[i, membership_cols, drop = TRUE])]
    paste(hits, collapse = " | ")
  }, character(1))

  out %>%
    dplyr::arrange(dplyr::desc(n_sets), gene)
}

# 1.6 Join Venn membership with DEG statistics from each comparison
build_venn_union_stats <- function(membership_df,
                                   detail_tables,
                                   feature_col = "feature") {
  out <- membership_df

  for (nm in names(detail_tables)) {
    df <- detail_tables[[nm]]
    if (is.null(df) || nrow(df) == 0) next

    df2 <- df %>%
      dplyr::rename(gene = dplyr::all_of(feature_col)) %>%
      dplyr::rename_with(~ paste0(nm, "__", .x), -gene)

    out <- out %>%
      dplyr::left_join(df2, by = "gene")
  }

  out
}

# 1.7 Plot and save a triple Venn diagram (Fig F/G style)
save_triple_venn <- function(sets_named_list,
                             title,
                             out_prefix,
                             fill = c("#e41a1c", "#a6cee3", "#a6cee3"),
                             alpha = 0.45) {

  venn_grob <- VennDiagram::venn.diagram(
    x = sets_named_list,
    filename = NULL,
    fill = fill,
    alpha = alpha,
    lwd = 2,
    col = "black",
    cex = 1.2,
    fontface = "bold",
    cat.cex = 1.1,
    cat.fontface = "bold",
    cat.dist = c(0.08, 0.08, 0.08),
    margin = 0.06,
    main = title,
    main.cex = 1.2,
    main.fontface = "bold"
  )

  ## PNG
  png(file.path(fig_dir, paste0(out_prefix, ".png")),
      width = 2400, height = 2000, res = 300)
  grid::grid.newpage()
  grid::grid.draw(venn_grob)
  dev.off()

  ## PDF
  pdf(file.path(fig_dir, paste0(out_prefix, ".pdf")),
      width = 8.5, height = 7)
  grid::grid.newpage()
  grid::grid.draw(venn_grob)
  dev.off()

  message("Saved Venn -> ", out_prefix)
}

# 1.8 Pick top N genes by |logFC| from a DEG table (for heatmap panels)
pick_top_by_abs_fc <- function(deg_df, genes, n = 8, feature_col = "feature") {
  if (length(genes) == 0) return(character(0))
  fc_col <- detect_fc_col(deg_df)

  deg_df %>%
    dplyr::filter(.data[[feature_col]] %in% genes) %>%
    dplyr::mutate(abs_fc = abs(.data[[fc_col]])) %>%
    dplyr::arrange(dplyr::desc(abs_fc), p_val_adj) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::pull(.data[[feature_col]]) %>%
    unique()
}

## Read DEG tables --------------------------------------------------------

deg_len <- read_step4_sheet(xlsx_len, sheet = "DEG_RNA")
deg_amg <- read_step4_sheet(xlsx_amg, sheet = "DEG_RNA")
deg_ko  <- read_step4_sheet(xlsx_ko,  sheet = "DEG_RNA")

## Build UP/DOWN gene sets per comparison
len_up <- get_deg_genes(deg_len, "up",   padj_cut, logfc_cut)
len_dn <- get_deg_genes(deg_len, "down", padj_cut, logfc_cut)

amg_up <- get_deg_genes(deg_amg, "up",   padj_cut, logfc_cut)
amg_dn <- get_deg_genes(deg_amg, "down", padj_cut, logfc_cut)

ko_up  <- get_deg_genes(deg_ko,  "up",   padj_cut, logfc_cut)
ko_dn  <- get_deg_genes(deg_ko,  "down", padj_cut, logfc_cut)

message("LEN_vs_CTR: Up=", length(len_up), " Down=", length(len_dn))
message("LENnAMG_vs_LEN: Up=", length(amg_up), " Down=", length(amg_dn))
message("LENnKO_vs_LEN: Up=", length(ko_up),  " Down=", length(ko_dn))

## Keep full DEG tables for export
fig7_len_up_detail <- get_deg_table(
  deg_len, "up", padj_cut, logfc_cut,
  comparison = "LEN_vs_CTR",
  venn_panel = "Fig7",
  set_label = "LEN_vs_CTR_Up"
)
fig7_amg_dn_detail <- get_deg_table(
  deg_amg, "down", padj_cut, logfc_cut,
  comparison = "LENnAMG_vs_LEN",
  venn_panel = "Fig7",
  set_label = "LENnAMG_vs_LEN_Dn"
)
fig7_ko_dn_detail <- get_deg_table(
  deg_ko, "down", padj_cut, logfc_cut,
  comparison = "LENnKO_vs_LEN",
  venn_panel = "Fig7",
  set_label = "LENnKO_vs_LEN_Dn"
)

fig8_len_dn_detail <- get_deg_table(
  deg_len, "down", padj_cut, logfc_cut,
  comparison = "LEN_vs_CTR",
  venn_panel = "Fig8",
  set_label = "LEN_vs_CTR_Dn"
)
fig8_amg_up_detail <- get_deg_table(
  deg_amg, "up", padj_cut, logfc_cut,
  comparison = "LENnAMG_vs_LEN",
  venn_panel = "Fig8",
  set_label = "LENnAMG_vs_LEN_Up"
)
fig8_ko_up_detail <- get_deg_table(
  deg_ko, "up", padj_cut, logfc_cut,
  comparison = "LENnKO_vs_LEN",
  venn_panel = "Fig8",
  set_label = "LENnKO_vs_LEN_Up"
)

## Venn diagrams ----------------------------------------------------------

## Fig 7: LEN_vs_CTR Up intersect with (LENnAMG_vs_LEN Down) and (LENnKO_vs_LEN Down)
venn_F <- list(
  "Len vs CTR\nUp"         = len_up,
  "Len&AMG vs Len\nDn"     = amg_dn,
  "Len&CKO vs Len\nDn"     = ko_dn
)
venn_F

save_triple_venn(
  sets_named_list = venn_F,
  title = "Up (LEN vs CTR) overlapped with Down (LENnAMG/LENnKO vs LEN)",
  out_prefix = "fig7_Venn_UP_reversal"
)


## Fig 8: LEN_vs_CTR Down intersect with (LENnAMG_vs_LEN Up) and (LENnKO_vs_LEN Up)
venn_G <- list(
  "Len vs CTR\nDn"         = len_dn,
  "Len&AMG vs Len\nUp"     = amg_up,
  "Len&CKO vs Len\nUp"     = ko_up
)

fig7_membership <- build_venn_membership(
  venn_F,
  membership_labels = c("LEN_vs_CTR_Up", "LENnAMG_vs_LEN_Dn", "LENnKO_vs_LEN_Dn")
)
fig8_membership <- build_venn_membership(
  venn_G,
  membership_labels = c("LEN_vs_CTR_Dn", "LENnAMG_vs_LEN_Up", "LENnKO_vs_LEN_Up")
)

save_triple_venn(
  sets_named_list = venn_G,
  title = "Down (LEN vs CTR) overlapped with Up (LENnAMG/LENnKO vs LEN)",
  out_prefix = "fig8_Venn_DN_reversal",
  fill = c("#377eb8", "#fdbf6f", "#fdbf6f")
)


## Export overlap gene lists ----------------------------------------------

## Overlap definitions (useful for “reversal” highlighting)
up_rev_amg  <- intersect(len_up, amg_dn)
up_rev_ko   <- intersect(len_up, ko_dn)
up_rev_both <- Reduce(intersect, list(len_up, amg_dn, ko_dn))

dn_rev_amg  <- intersect(len_dn, amg_up)
dn_rev_ko   <- intersect(len_dn, ko_up)
dn_rev_both <- Reduce(intersect, list(len_dn, amg_up, ko_up))

overlap_sheets <- list(
  "LEN_up"             = data.frame(gene = sort(len_up)),
  "LEN_down"           = data.frame(gene = sort(len_dn)),
  "AMG_down_for_Fig7"  = data.frame(gene = sort(amg_dn)),
  "KO_down_for_Fig7"   = data.frame(gene = sort(ko_dn)),
  "AMG_up_for_Fig8"    = data.frame(gene = sort(amg_up)),
  "KO_up_for_Fig8"     = data.frame(gene = sort(ko_up)),
  "Up_reversed_by_AMG" = data.frame(gene = sort(up_rev_amg)),
  "Up_reversed_by_KO"  = data.frame(gene = sort(up_rev_ko)),
  "Up_reversed_by_BOTH"= data.frame(gene = sort(up_rev_both)),
  "Down_reversed_by_AMG"=data.frame(gene = sort(dn_rev_amg)),
  "Down_reversed_by_KO" =data.frame(gene = sort(dn_rev_ko)),
  "Down_reversed_by_BOTH"=data.frame(gene = sort(dn_rev_both))
)

out_overlap_xlsx <- file.path(out_dir, "Step5_Venn_overlap_gene_lists.xlsx")
openxlsx2::write_xlsx(overlap_sheets, file = out_overlap_xlsx)
message("Saved overlap gene lists -> ", out_overlap_xlsx)

fig7_union_stats <- build_venn_union_stats(
  membership_df = fig7_membership,
  detail_tables = list(
    LEN_vs_CTR_Up = fig7_len_up_detail,
    LENnAMG_vs_LEN_Dn = fig7_amg_dn_detail,
    LENnKO_vs_LEN_Dn = fig7_ko_dn_detail
  )
)

fig8_union_stats <- build_venn_union_stats(
  membership_df = fig8_membership,
  detail_tables = list(
    LEN_vs_CTR_Dn = fig8_len_dn_detail,
    LENnAMG_vs_LEN_Up = fig8_amg_up_detail,
    LENnKO_vs_LEN_Up = fig8_ko_up_detail
  )
)

venn_deg_info <- tibble::tibble(
  padj_cutoff = padj_cut,
  logfc_cutoff = logfc_cut,
  len_source_xlsx = xlsx_len,
  amg_source_xlsx = xlsx_amg,
  ko_source_xlsx = xlsx_ko,
  fig7_output = "fig7_Venn_UP_reversal",
  fig8_output = "fig8_Venn_DN_reversal"
)

venn_deg_sheets <- list(
  "Info" = venn_deg_info,
  "Fig7_LEN_Up" = fig7_len_up_detail,
  "Fig7_AMG_Dn" = fig7_amg_dn_detail,
  "Fig7_KO_Dn" = fig7_ko_dn_detail,
  "Fig7_Membership" = fig7_membership,
  "Fig7_UnionStats" = fig7_union_stats,
  "Fig8_LEN_Dn" = fig8_len_dn_detail,
  "Fig8_AMG_Up" = fig8_amg_up_detail,
  "Fig8_KO_Up" = fig8_ko_up_detail,
  "Fig8_Membership" = fig8_membership,
  "Fig8_UnionStats" = fig8_union_stats
)

out_venn_deg_xlsx <- file.path(out_dir, "Step5_Venn_DEG_details.xlsx")
openxlsx2::write_xlsx(venn_deg_sheets, file = out_venn_deg_xlsx)
message("Saved Venn DEG details -> ", out_venn_deg_xlsx)





### Section 8: Figure 9. Key-gene heatmap -------------------------------------
## Heatmap for key genes (SCT assay, group-averaged, z-scored)
## Input : Seurat object "seu"
## Output: figures/final/Fig9_SCT_keygenes_heatmap.png/pdf

seu <- seu_full[,seu_full$group %in% c("CTR", "LEN","LENnAMG", "LENnKO")]
seu <- seu[, seu$CellType_l1 == "Neu"]


## -------------------------- Parameters -----------------------------------
out_png <- file.path(workdir, "figures", "final", "Fig9_SCT_keygenes_heatmap.png")
out_pdf <- file.path(workdir, "figures", "final", "Fig9_SCT_keygenes_heatmap.pdf")
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

# group order (paper-like: Len, Len&CKO, Len&AMG)
group_order <- c("CTR" ,"LEN","LENnAMG", "LENnKO")

# heatmap color scale (same style you used)
hm_colors5 <- c("#1f618d", "#2E86C1", "white", "#ec7063", "#AD272B")

## -------------------------- Gene modules ---------------------------------
gene_modules <- list(
  `PMN-MDSC signatures` = c(
    "S100a8","S100a9","Lcn2","Camp","Retnlg",
    "Ngp","Ly6g","Ly6c2","Lyz2"
  ),
  `Myeloid activation` = c(
    "Padi4","Mmp8","Pglyrp1","Serpinb1a","Hp","Cebpd",
    "Tnfaip8l2","Cd55","Ceacam1","Cd177"
  ),
  `Chemotaxis / Extravasation` = c(
    "Alox5ap","Ltb4r1","Selp",
    "Tagln2","Capg","Flna","Arpc3",
    "Dock8","Rasgrp2","Rgs18","Itgb2l"
  ),
  `IFN-stimulated genes` = c(
    "Stat1","Irf1","Jak2","Tbk1",
    "Nfkb1","Irak2",
    "Cxcl10","Gbp2","Gbp3","Gbp5",
    "Ifi47","Ifi207","Igtp","Irgm1","Ly6e"
  ),
  `Antigen processing & presentation` = c(
    "B2m","H2-K1",
    "Tap1","Tap2","Tapbp",
    "Nlrc5","Psme2","Psme2b",
    "Ctss","Cd86"
  )
)

gene_order <- unlist(gene_modules, use.names = FALSE)
gene_module <- rep(names(gene_modules), times = lengths(gene_modules))
names(gene_module) <- gene_order

## -------------------------- Checks + subset ------------------------------
stopifnot("group" %in% colnames(seu@meta.data))
stopifnot("SCT" %in% names(seu@assays))

groups_use <- intersect(group_order, unique(as.character(seu$group)))
if (length(groups_use) < 2) stop("Not enough groups found in seu$group. Found: ", paste(groups_use, collapse = ", "))

seu_sub <- subset(seu, subset = group %in% groups_use)
seu_sub$group <- factor(as.character(seu_sub$group), levels = groups_use)

Seurat::DefaultAssay(seu_sub) <- "SCT"

## -------------------------- Average expression (SCT) ----------------------
avg <- Seurat::AverageExpression(
  seu_sub,
  assays   = "SCT",
  slot     = "data",
  group.by = "group",
  verbose  = FALSE
)

mat <- avg$SCT  # genes x groups

# keep gene order as you provided
gene_order_use <- gene_order[gene_order %in% rownames(mat)]
missing_genes  <- setdiff(gene_order, gene_order_use)
if (length(missing_genes) > 0) {
  message("Missing genes in SCT assay (will be dropped): ", paste(missing_genes, collapse = ", "))
}

mat_sub <- mat[gene_order_use, groups_use, drop = FALSE]  # genes x groups

# z-score per gene across groups + clamp [-2, 2]
mat_z <- t(scale(t(mat_sub)))
mat_z[is.na(mat_z)] <- 0
mat_z <- pmax(pmin(mat_z, 2), -2)

mat_plot <- t(mat_z)  # groups x genes (rows=groups, cols=genes)

# module factor for columns
module_factor <- factor(gene_module[gene_order_use], levels = names(gene_modules))

## -------------------------- Reorder genes within each module --------------
# We reorder genes to improve visual continuity.
# Criterion: delta = Z(LEN) - Z(CTR). Larger delta => more LEN-high.
# For IFN/Antigen modules (often LEN-low), we sort ascending to put LEN-low genes first.

stopifnot(all(c("CTR","LEN") %in% rownames(mat_plot)))  # needed for delta

delta_len_ctr <- mat_plot["LEN", ] - mat_plot["CTR", ]  # named vector (genes)

# choose per-module sorting direction
module_sort_dir <- c(
  `PMN-MDSC signatures`               = "desc",  # LEN-high first
  `Myeloid activation`                = "desc",
  `Chemotaxis / Extravasation`        = "desc",
  `IFN-stimulated genes`              = "asc",   # LEN-low first
  `Antigen processing & presentation` = "asc"    # LEN-low first
)

# build a reordered gene list
gene_order_reordered <- unlist(lapply(names(gene_modules), function(mod) {
  genes_mod <- intersect(gene_modules[[mod]], colnames(mat_plot))
  if (length(genes_mod) == 0) return(character(0))

  d <- delta_len_ctr[genes_mod]
  d <- d[!is.na(d)]

  if (module_sort_dir[[mod]] == "asc") {
    names(sort(d, decreasing = FALSE))
  } else {
    names(sort(d, decreasing = TRUE))
  }
}), use.names = FALSE)

# apply order to matrix + module_factor
mat_plot <- mat_plot[, gene_order_reordered, drop = FALSE]
module_factor <- factor(gene_module[gene_order_reordered], levels = names(gene_modules))

## -------------------------- Plot (ComplexHeatmap preferred) ---------------
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# heat colors
col_fun <- circlize::colorRamp2(
  c(-2, -1, 0, 1, 2),
  hm_colors5
)

# Module block colors for the top annotation.
module_colors <- c(
  `PMN-MDSC signatures` = "#f4cccc",
  `Myeloid activation` = "#fce5cd",
  `Chemotaxis / Extravasation` = "#d9ead3",
  `IFN-stimulated genes` = "#cfe2f3",
  `Antigen processing & presentation` = "#ead1dc"
)

top_anno <- ComplexHeatmap::HeatmapAnnotation(
  Module = ComplexHeatmap::anno_block(
    gp = grid::gpar(fill = module_colors[names(gene_modules)], col = NA),
    labels = names(gene_modules),
    labels_gp = grid::gpar(fontface = "bold", fontsize = 11)
  )
)

ht <- ComplexHeatmap::Heatmap(
  mat_plot,
  name = "Z",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  column_names_rot = 45,
  column_split = module_factor,
  top_annotation = top_anno,
  heatmap_legend_param = list(title = "Z", at = c(-2, -1, 0, 1, 2))
)

# PNG
png(out_png, width = 4200, height = 900, res = 300)
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
dev.off()

# PDF
pdf(out_pdf, width = 18, height = 4.2)
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
dev.off()

message("Saved: ", out_png)
message("Saved: ", out_pdf)
