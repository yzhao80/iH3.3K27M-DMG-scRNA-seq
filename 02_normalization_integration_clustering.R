##=============================================================================
## 02_normalization_integration_clustering.R
##
## Project : iH3.3K27M DMG scRNA-seq (Khairkhah et al.)
## Purpose : Log-normalization, feature selection, PCA, RPCA batch-correction,
##           graph-based clustering, and UMAP embedding.
## Methods : Seurat v5.1.0 – LogNormalize, VST, RPCA integration, Louvain
##           clustering (res = 0.4), RunUMAP (12 dims)
## Input   : seurat_merged_QC.rds
## Output  : seurat_merged_clustered.rds  (annotated clusters, UMAP)
##=============================================================================

## ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(clustree)      # resolution stability
library(RColorBrewer)

## ── Paths ─────────────────────────────────────────────────────────────────────
input_dir  <- "results/01_QC/"
output_dir <- "results/02_clustering/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ── 1. Load QC-passed object ─────────────────────────────────────────────────
seurat_merged <- readRDS(file.path(input_dir, "seurat_merged_QC.rds"))

## ── 2. Normalization, variable features, scaling, PCA ────────────────────────
seurat_merged <- NormalizeData(seurat_merged,
    normalization.method = "LogNormalize", scale.factor = 1e4)
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst")
seurat_merged <- ScaleData(seurat_merged)   # scale all genes for DE later
seurat_merged <- RunPCA(seurat_merged,
    features = VariableFeatures(seurat_merged), npcs = 50)

# Elbow plot to determine dimensionality – 12 PCs selected
p_elbow <- ElbowPlot(seurat_merged, ndims = 50) +
    geom_vline(xintercept = 12, color = "red", linetype = "dashed") +
    labs(title = "Elbow plot – 12 PCs selected")
ggsave(file.path(output_dir, "elbow_plot.pdf"), p_elbow, width = 6, height = 4)

## ── 3. RPCA integration across four samples ───────────────────────────────────
# Split layers by sample, then integrate to correct batch effects
seurat_merged[["RNA"]] <- split(seurat_merged[["RNA"]],
                                f = seurat_merged$sample_id)

seurat_merged <- IntegrateLayers(
    object        = seurat_merged,
    method        = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction  = "integrated.rpca",
    verbose        = TRUE
)

## ── 4. Neighbor graph, clustering, and UMAP ──────────────────────────────────
# All downstream neighbour-finding and UMAP use the RPCA-integrated embedding
N_DIMS <- 12

seurat_merged <- FindNeighbors(seurat_merged,
    reduction = "integrated.rpca", dims = 1:N_DIMS)

# Cluster across a range of resolutions; res = 0.4 yielded 25 biologically
# interpretable clusters
seurat_merged <- FindClusters(seurat_merged,
    resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0))

seurat_merged <- RunUMAP(seurat_merged,
    reduction = "integrated.rpca", dims = 1:N_DIMS)

# Resolution stability (clustree)
res_cols <- grep("_snn_res\\.", colnames(seurat_merged@meta.data), value = TRUE)
pref     <- sub("(.*_snn_res\\.).*", "\\1", res_cols[1])
p_tree   <- clustree(seurat_merged@meta.data, prefix = pref)
ggsave(file.path(output_dir, "clustree_resolution_stability.pdf"),
       p_tree, width = 8, height = 10)

# Fix active identity to res = 0.4
Idents(seurat_merged) <- "RNA_snn_res.0.4"

## ── 5. UMAP visualizations ───────────────────────────────────────────────────
p_umap_clusters <- DimPlot(seurat_merged, reduction = "umap",
    label = TRUE, repel = TRUE) + NoLegend()
p_umap_treatment <- DimPlot(seurat_merged,
    group.by = "treatment",
    split.by = "treatment", ncol = 2, label = FALSE)
p_umap_sample <- DimPlot(seurat_merged,
    group.by = "sample_id",
    split.by = "sample_id", ncol = 2, label = FALSE)

ggsave(file.path(output_dir, "UMAP_clusters.pdf"),
       p_umap_clusters, width = 7, height = 6)
ggsave(file.path(output_dir, "UMAP_by_treatment.pdf"),
       p_umap_treatment, width = 12, height = 6)
ggsave(file.path(output_dir, "UMAP_by_sample.pdf"),
       p_umap_sample, width = 12, height = 6)

## ── 6. Cell-type annotation ───────────────────────────────────────────────────
# Manual annotation based on canonical marker expression (see 03_DE.R for
# marker identification). Final labels used in all downstream analyses.
cluster_labels <- c(
    "0"  = "Neurons",
    "1"  = "Neurons",
    "2"  = "Neurons",
    "3"  = "Excitatory Neurons",
    "4"  = "Neurons",
    "5"  = "Neurons",
    "6"  = "Neurons",
    "7"  = "Inhibitory Neurons",
    "8"  = "Neurons",
    "9"  = "Astrocyte",
    "10" = "Microglia",
    "11" = "Endothelial Cells",
    "12" = "Fibroblasts",
    "13" = "Ependymal Cells",
    "14" = "Epithelial Cells",
    "15" = "Immune Cells",
    "16" = "DMG Tumor",   # PDGFRA+ OPC (Early)
    "17" = "DMG Tumor",   # PDGFRA+ OPC (Differentiating)
    "18" = "Purkinje Neurons",
    "19" = "Oligodendrocytes",
    "20" = "Inhibitory Neurons",
    "21" = "Neurons",
    "22" = "Neurons",
    "23" = "Excitatory Neurons",
    "24" = "Neurons"
)

seurat_merged <- RenameIdents(seurat_merged, cluster_labels)
seurat_merged$CellType_alpha <- factor(
    as.character(Idents(seurat_merged)),
    levels = sort(unique(as.character(Idents(seurat_merged))))
)

cat("\nCell-type composition:\n")
print(table(seurat_merged$CellType_alpha, seurat_merged$treatment))

## ── 7. Annotated UMAP (Fig. 5E) ──────────────────────────────────────────────
# Figure 5E: UMAP colored by cell type, split by ON/OFF condition
mycol <- c(
    "Astrocyte"          = "#AEC6CF",
    "DMG Tumor"          = "#1F78B4",
    "Endothelial Cells"  = "#B2DF8A",
    "Ependymal Cells"    = "#33A02C",
    "Epithelial Cells"   = "#FB9A99",
    "Excitatory Neurons" = "#E31A1C",
    "Fibroblasts"        = "#FDBF6F",
    "Immune Cells"       = "#FF7F00",
    "Inhibitory Neurons" = "#CAB2D6",
    "Microglia"          = "#6A3D9A",
    "Neurons"            = "#FFFF99",
    "Oligodendrocytes"   = "#B15928",
    "Purkinje Neurons"   = "#A6CEE3"
)

seurat_merged@meta.data$treatment <- factor(
    seurat_merged@meta.data$treatment, levels = c("ON", "OFF"))

p_annotated <- DimPlot(
    seurat_merged,
    reduction  = "umap",
    group.by   = "CellType_alpha",
    split.by   = "treatment",
    cols       = mycol,
    label      = TRUE,
    label.size = 3,
    repel      = TRUE,
    ncol       = 2
) + theme(legend.position = "right")

ggsave(file.path(output_dir, "UMAP_annotated_ON_OFF_Fig5E.pdf"),
       p_annotated, width = 14, height = 6)

## ── 8. Cell-type proportion bar chart (Fig. 5F) ───────────────────────────────
# Figure 5F: Stacked bar chart showing cell-type proportions per condition
prop_df <- as.data.frame(
    prop.table(table(seurat_merged$CellType_alpha, seurat_merged$treatment),
               margin = 2)
)
colnames(prop_df) <- c("CellType", "Treatment", "Proportion")

p_bar <- ggplot(prop_df, aes(x = Treatment, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = mycol) +
    labs(x = NULL, y = "Cell-type proportion", fill = "CellType") +
    theme_classic(base_size = 12)

ggsave(file.path(output_dir, "CellType_proportion_bar_Fig5F.pdf"),
       p_bar, width = 5, height = 6)

## ── 9. Save ───────────────────────────────────────────────────────────────────
saveRDS(seurat_merged, file.path(output_dir, "seurat_merged_clustered.rds"))
cat("Saved:", file.path(output_dir, "seurat_merged_clustered.rds"), "\n")
