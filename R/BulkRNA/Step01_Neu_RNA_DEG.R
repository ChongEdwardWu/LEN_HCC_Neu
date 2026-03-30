### Section 0: Preparation --------------------------------------------------

# clear the environment
rm(list = ls())
gc()

# Optional: load a personal R profile if you use one locally
if (file.exists("path_to_.radian_profile")) source("path_to_.radian_profile")
# .libPaths()

# set your working directory
workdir <- "path_to_bulk_rna_analysis"
setwd(workdir)
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

library("tximeta")
library(tidyverse)
library(openxlsx2)
library(AnnotationDbi)
library(org.Mm.eg.db)

# ! httpgd::hgd()

### Section 1: Summarize transcript-level quantifications to gene-level --------------------------------------------------

# Data import
gse_rds <- file.path(workdir, "HCC_trt_neu.rds")
stopifnot(file.exists(gse_rds))
gse <- readRDS(gse_rds)

gse@colData$position[7] <- "TU"
gse.TU <- subset(gse,,(gse$treatment == "Lenva" | gse$treatment == "Untx" )& gse$position == "TU")
gse.TU <- subset(gse.TU,,(gse.TU$names != "Untx-TU-PMN-2"))
gse.TU$animal <- c("1","2","3","4")






### Section 2: Gene Expression Normalization to TPM --------------------------------------------------
# Extract the TPM (abundance) matrix from the gene-level SummarizedExperiment (gse.TU)
TPM <- gse.TU@assays@data$abundance
# row names should be shortened before mapping
ens_TPM <- substr(rownames(TPM), 1, 18)
sym_TPM <- mapIds(org.Mm.eg.db,
                  keys = ens_TPM,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")
TPM_df <- TPM %>%
  as_tibble %>%
  mutate(GeneSymbol = sym_TPM) %>%
  mutate(Ensembl = ens_TPM) %>%
  filter(!is.na(GeneSymbol)) %>%
  dplyr::select(Ensembl, GeneSymbol, everything()) %>%
  arrange(Ensembl)

### Section 3: Export Gene Expression Counts --------------------------------------------------
# Extract the raw counts matrix from gse.TU
counts <- gse.TU@assays@data$counts
# row names should be shortened before mapping
ens_counts <- substr(rownames(counts), 1, 18)
sym_counts <- mapIds(org.Mm.eg.db,
                     keys = ens_counts,
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")
counts_df <- counts %>%
  as_tibble %>%
  mutate(GeneSymbol = sym_counts) %>%
  mutate(Ensembl = ens_counts) %>%
  dplyr::filter(!is.na(GeneSymbol)) %>%
  dplyr::select(Ensembl, GeneSymbol, everything()) %>%
  arrange(Ensembl)

# Export both TPM and Counts into the same Excel file with two sheets
resultfile <- file.path("results", "NeuBulk_LEN_vs_CON_results.xlsx")
wb <- openxlsx2::wb_workbook()
wb$add_worksheet("counts")
wb$add_data("counts", x = as.data.frame(counts_df))
wb$add_worksheet("TPM")
wb$add_data("TPM", x = as.data.frame(TPM_df))
wb$save(resultfile, overwrite = TRUE)

### Section 4: DESeqDataSet Object Construction and Exploratory Analysis --------------------------------------------------
library(DESeq2)

gse.TU$treatment <- droplevels(as.factor(gse.TU$treatment))
table(gse.TU$treatment)

# Ensure that the colData of gse.TU contains a column named "treatment" with values "CON" and "LEN"
dds <- DESeqDataSet(gse.TU, design = ~ treatment)
# Pre-filter lowly expressed genes (keep genes with a total count > 5)
keep <- rowSums(counts(dds)) > 5
dds <- dds[keep, ]



### Section 5: Differential Expression Analysis Using DESeq2 --------------------------------------------------
# Run the DESeq2 pipeline
dds <- DESeq(dds)
# Compare LEN vs CON (log2 fold change of LEN relative to CON)
res <- results(dds, contrast = c("treatment", "Lenva", "Untx"), alpha = 0.05)
summary(res)

# Annotate the DE results: map Ensembl IDs to gene symbols and Entrez IDs
# annotate results
library("AnnotationDbi")
library("org.Mm.eg.db")
# row names should be shortened before mapping
ens_ids <- substr(rownames(res), 1, 18)
res$GeneSymbol <- mapIds(org.Mm.eg.db,
                         keys = ens_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
res$Entrez <- mapIds(org.Mm.eg.db,
                     keys = ens_ids,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")
# Order the results by adjusted p-value
res_ord <- as.data.frame(res[order(res$padj), ]) %>%
  dplyr::filter(!is.na(GeneSymbol))

build_bad_features <- function(gene_names) {
  hist_genes    <- grep("^Hist", gene_names, ignore.case = TRUE, value = TRUE)
  hb_genes      <- grep("^Hb[ab]-|^HB[^(P)]", gene_names, value = TRUE)
  mt_genes      <- grep("^mt-|^MT-|MTRNR2L|Mtrnr2l|^Mtmr", gene_names, ignore.case = TRUE, value = TRUE)
  rps_genes     <- grep("^Rp[sl]|^RP[SL]", gene_names, ignore.case = TRUE, value = TRUE)
  rik_genes     <- grep("^.*Rik$|^Rik", gene_names, ignore.case = TRUE, value = TRUE)
  pseudo_genes  <- grep("-rs|-ps", gene_names, ignore.case = TRUE, value = TRUE)
  mir_genes     <- grep("^Mir", gene_names, ignore.case = TRUE, value = TRUE)
  gencode_genes <- grep("^Gm|ENSMUSG", gene_names, ignore.case = TRUE, value = TRUE)

  unique(c(hist_genes, hb_genes, mt_genes, rps_genes, rik_genes, pseudo_genes, mir_genes, gencode_genes))
}

bad_features <- build_bad_features(res_ord$GeneSymbol)
good_genes <- setdiff(res_ord$GeneSymbol, bad_features)

# filter genes
res_ord <- res_ord[res_ord$GeneSymbol %in% good_genes,]

# Export the DE results to the xlsx file
wb$add_worksheet("DEG")
wb$add_data("DEG", x = as.data.frame(res_ord))
wb$save(resultfile, overwrite = TRUE)

### Section 6: GO Over-representation Analysis ----------------------
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(tidyverse)

deg <- as.data.frame(res_ord)

# Filter for upregulated genes: log2FoldChange > 1, padj < 0.05 and non-missing Entrez IDs
up_genes <- deg %>%
  filter(log2FoldChange > 1, padj < 0.05, !is.na(Entrez)) %>%
  pull(Entrez) %>%
  as.character()

# Define the universe as all genes in the DEG results with non-missing Entrez IDs
universe <- deg %>%
  filter(!is.na(Entrez)) %>%
  pull(Entrez) %>%
  as.character()

# Perform GO enrichment analysis for Biological Process (BP) using enrichGO
egoup <- enrichGO(gene          = up_genes,
                  universe      = universe,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  minGSSize     = 50,
                  maxGSSize     = 200,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

# Improve visualization by computing pairwise term similarity
egoupfilt <- pairwise_termsim(egoup)
# Export the filtered up-regulated GO enrichment results to the result Excel file
wb$add_worksheet("GO_BP_upfilt")
wb$add_data("GO_BP_upfilt", x = as.data.frame(egoupfilt@result))
wb$save(resultfile, overwrite = TRUE)

# # Create a dotplot of the enriched GO terms
# dotplot_up_GO <- dotplot(egoupfilt, showCategory = 10) +
#   ggtitle("GO Enrichment Analysis (Biological Process)")
# ggsave(filename = "figures/Step01_Bulk_LEN_vs_CON_GO_up_dotplot.pdf",
#        plot = dotplot_up_GO,
#        width = 7,
#        height = 7,
#        units = "in")


# Filter for upregulated genes: log2FoldChange < -1, padj < 0.05 and non-missing Entrez IDs
dn_genes <- deg %>%
  filter(log2FoldChange < -1, padj < 0.05, !is.na(Entrez)) %>%
  pull(Entrez) %>%
  as.character()

# Define the universe as all genes in the DEG results with non-missing Entrez IDs
universe <- deg %>%
  filter(!is.na(Entrez)) %>%
  pull(Entrez) %>%
  as.character()

# # Perform GO enrichment analysis for Biological Process (BP) using enrichGO
egodn <- enrichGO(gene          = dn_genes,
                  universe      = universe,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  minGSSize     = 50,
                  maxGSSize     = 200,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

# Improve visualization by computing pairwise term similarity
egodnfilt <- pairwise_termsim(egodn)
# Export the filtered down-regulated GO enrichment results to the result Excel file
wb$add_worksheet("GO_BP_dnfilt")
wb$add_data("GO_BP_dnfilt", x = as.data.frame(egodnfilt@result))
wb$save(resultfile, overwrite = TRUE)

# Create a dotplot of the enriched GO terms
# dotplot_dn_GO <- dotplot(egodnfilt, showCategory = 15) +
#   ggtitle("GO Enrichment Analysis (Biological Process)")
# ggsave(filename = "figures/Step01_Bulk_LEN_vs_CON_GO_dn_dotplot.pdf",
#        plot = dotplot_dn_GO,
#        width = 7,
#        height = 7,
#        units = "in")


### Section 7: Gene Set Enrichment Analysis for Bulk LEN vs CON ------------------------------------------------------------
## Load packages and prepare environment
library(clusterProfiler)
library(org.Mm.eg.db)    # for mouse
library(stringr)
library(fgsea)
library(msigdbr)
set.seed(123)
library(future)
nworkers <- 16

# List available gene set categories for mouse (for reference)
m_df <- msigdbr(species = "Mus musculus")
print(m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat), n = 23)

## Prepare gene list for gsea
# Read DEG results from the Excel file (DEG sheet)
# resultfile <- "results/Bulk_LEN_vs_CON_results.xlsx"
DEG <- openxlsx2::read_xlsx(resultfile, sheet = "DEG", row.names = 1)
# View(DEG)


# Extract log2 fold change values and assign gene names from the "GeneSymbol" column
logFC <- DEG$log2FoldChange
names(logFC) <- DEG$GeneSymbol
# Sorting gene list in decreasing order of logFC
geneList <- sort(logFC, decreasing = TRUE)
# Remove genes without names
geneList <- geneList[!is.na(names(geneList))]
# Display the length of the gene list
length(geneList)

## gsea for Reactome Pathways (C2:CP:REACTOME)
pathwaysDF <- msigdbr("mouse", category = "C2", subcategory = "CP:REACTOME")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)
length(pathways)
fgseaRes <- fgsea(pathways, geneList, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
REACTOMEres <- fgseaRes
head(REACTOMEres)
wb$add_worksheet("gsea_REACTOME")
wb$add_data("gsea_REACTOME", x = as.data.frame(REACTOMEres))
wb$save(resultfile, overwrite = TRUE)

## gsea for KEGG Pathways (C2:CP:KEGG)
pathwaysDF <- msigdbr("mouse", category = "C2", subcategory = "CP:KEGG")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)
length(pathways)
fgseaRes <- fgsea(pathways, geneList, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
KEGGres <- fgseaRes
head(KEGGres)
wb$add_worksheet("gsea_KEGG")
wb$add_data("gsea_KEGG", x = as.data.frame(KEGGres))
wb$save(resultfile, overwrite = TRUE)



## gsea for GO Biological Process (C5:GO:BP)
pathwaysDF <- msigdbr("mouse", category = "C5", subcategory = "GO:BP")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)
fgseaRes <- fgsea(pathways, geneList, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
GOBPres <- fgseaRes
head(GOBPres)

wb$add_worksheet("gsea_GO_BP")
wb$add_data("gsea_GO_BP", x = as.data.frame(GOBPres))
wb$save(resultfile, overwrite = TRUE)

## gsea for Hallmark Gene Sets
pathwaysDF <- msigdbr("mouse", category = "H")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)
fgseaRes <- fgsea(pathways, geneList, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
Hres <- fgseaRes
head(Hres)

wb$add_worksheet("gsea_HallMark")
wb$add_data("gsea_HallMark", x = as.data.frame(Hres))
wb$save(resultfile, overwrite = TRUE)


## gsea for Neutrophils-related Gene Sets
pathways <- readRDS("path_to_neutrophil_gene_sets_rds")
fgseaRes <- fgsea(pathways, geneList, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
Neures <- fgseaRes
head(Neures)
wb$add_worksheet("gsea_Neu")
wb$add_data("gsea_Neu", x = as.data.frame(Neures))
wb$save(resultfile, overwrite = TRUE)


# Figure 1: GO over-representation analysis dotplots
library(patchwork)
library(ggplot2)

## Create a dotplot of the enriched GO terms
dotplot_up_GO <- dotplot(egoupfilt, showCategory = 10) +
  ggtitle("Up-regulation")

# Create a dotplot of the enriched GO terms
dotplot_dn_GO <- dotplot(egodnfilt, showCategory = 10) +
  ggtitle("Down-regulation")


# Arrange the two panels horizontally.
p_GO_up_down <- (dotplot_up_GO | dotplot_dn_GO) 

# Optionally add an overall title.
p_GO_up_down <- p_GO_up_down +
  plot_annotation(title = "GO Over-representation Analysis (BP): LEN vs CON")

# Display the combined plot.
p_GO_up_down

# Save the combined plot.
ggsave("figures/Fig1_GO_BP_ORA_UP_DOWN_dotplot.png", p_GO_up_down, width = 350, height = 250, units = "mm", dpi = 300)
ggsave("figures/Fig1_GO_BP_ORA_UP_DOWN_dotplot.pdf", p_GO_up_down, width = 350, height = 250, units = "mm")

# Optional: Save session information for reproducibility
writeLines(capture.output(sessionInfo()), "session_info.txt")
