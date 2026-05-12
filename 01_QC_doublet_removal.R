##=============================================================================
## 01_QC_doublet_removal.R
##
## Project : iH3.3K27M DMG scRNA-seq (Khairkhah et al.)
## Purpose : Ambient RNA decontamination, doublet detection, and quality
##           filtering of the merged raw Seurat object produced by Cell Ranger.
## Methods : DecontX (celda v1.20.0), scDblFinder (v1.18.0), Seurat v5.1.0
## Input   : seurat_merged_raw.rds  – merged raw feature-barcode matrix
## Output  : seurat_merged_QC.rds   – filtered, QC-passed Seurat object
##=============================================================================

## ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)       # v5.1.0
library(celda)        # v1.20.0  – DecontX ambient-RNA removal
library(scDblFinder)  # v1.18.0  – doublet detection
library(ggplot2)
library(cowplot)
library(dplyr)

## ── Paths (edit to match your environment) ───────────────────────────────────
input_dir  <- "data/"
output_dir <- "results/01_QC/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ── 1. Load merged raw Seurat object ─────────────────────────────────────────
seurat_merged <- readRDS(file.path(input_dir, "seurat_merged_raw.rds"))

# Sample and treatment metadata
# Samples: NK13, NK14 = H3.3K27M ON; NK15, NK16 = H3.3K27M OFF
sample_ids <- c("NK13", "NK14", "NK15", "NK16")
treatment  <- c("ON",   "ON",   "OFF",  "OFF")
sex        <- c("F",    "M",    "F",    "M")

sample_map  <- setNames(treatment, sample_ids)
sample_map2 <- setNames(sex,       sample_ids)

seurat_merged@meta.data$treatment <- sample_map[seurat_merged@meta.data$sample_id]
seurat_merged@meta.data$sex       <- sample_map2[seurat_merged@meta.data$sample_id]

## ── 2. Mitochondrial gene content ────────────────────────────────────────────
seurat_merged[["percent.mt"]] <- PercentageFeatureSet(
    seurat_merged, pattern = "^(MT|mt)-"
)

## ── 3. DecontX – ambient RNA removal ─────────────────────────────────────────
# Run per sample to avoid batch confounding
seurat_merged <- decontX(seurat_merged)

## ── 4. scDblFinder – doublet detection ───────────────────────────────────────
set.seed(42)
sce <- as.SingleCellExperiment(seurat_merged)
sce <- scDblFinder(sce, samples = "sample_id")
seurat_merged$scDblFinder.class <- sce$scDblFinder.class
seurat_merged$scDblFinder.score <- sce$scDblFinder.score

# Remove doublets
seurat_merged <- subset(seurat_merged, subset = scDblFinder.class == "singlet")

## ── 5. QC plots (pre-filter) ─────────────────────────────────────────────────
p1 <- VlnPlot(seurat_merged, features = "nFeature_RNA", pt.size = 0) +
    geom_hline(yintercept = c(500, 10000), color = "red") +
    NoLegend()
p2 <- VlnPlot(seurat_merged, features = "nCount_RNA", pt.size = 0) +
    geom_hline(yintercept = 20000, color = "red") +
    NoLegend()
p3 <- VlnPlot(seurat_merged, features = "percent.mt", pt.size = 0) +
    geom_hline(yintercept = 5, color = "red") +
    NoLegend()

p_QC <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
ggsave(file.path(output_dir, "QC_violin_prefilter.pdf"),
       plot = p_QC, width = 14, height = 5)

## ── 6. Quality thresholding ───────────────────────────────────────────────────
# Thresholds: nFeature 500–10,000 | nCount < 20,000 | percent.mt < 5 %
seurat_merged <- subset(
    seurat_merged,
    subset = nFeature_RNA > 500 &
             nFeature_RNA < 10000 &
             nCount_RNA   < 20000 &
             percent.mt   < 5
)

cat("Cells after QC filtering:", ncol(seurat_merged), "\n")

## ── 7. Save ───────────────────────────────────────────────────────────────────
saveRDS(seurat_merged, file.path(output_dir, "seurat_merged_QC.rds"))
cat("Saved:", file.path(output_dir, "seurat_merged_QC.rds"), "\n")
