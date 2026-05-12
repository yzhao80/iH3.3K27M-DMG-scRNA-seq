##=============================================================================
## 07_immune_subclustering.R
##
## Project : iH3.3K27M DMG scRNA-seq (Khairkhah et al.)
## Purpose : Subset and re-cluster the "Immune Cells" population identified
##           in script 02, annotate immune subtypes, and generate the
##           publication-ready UMAP and composition barplot for ON vs OFF.
## Methods : Seurat v5.1.0. Immune cells (cluster 18 at res = 0.4) are
##           re-normalized, re-clustered independently (PCA 15 dims,
##           Louvain res = 0.3) and annotated using canonical immune markers.
##           Subtype composition compared ON vs OFF using Fisher's exact test.
## Input   : NK_seurat_merged.rds  (annotated object from script 02)
## Output  : NK_immune.rds                        – annotated immune Seurat object
##           immune_UMAP_ON_OFF_Fig6D.png          – Fig 6D
##           immune_composition_barplot_Fig6E.png  – Fig 6E
##           markers_immune_subclusters.csv        – subcluster marker genes
##           deg_ON_OFF_immune_all.csv             – ON vs OFF DEGs, all immune
##           deg_ON_OFF_immune_subcluster1.csv     – ON vs OFF DEGs, subcluster 1
## Figures : Fig 6D (immune UMAP, split by ON/OFF)
##           Fig 6E (immune cell-type composition barplot)
##=============================================================================

## ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(cowplot)
library(patchwork)
library(openxlsx)
library(scales)

## ── Paths ─────────────────────────────────────────────────────────────────────
input_dir  <- "results/02_clustering/"
output_dir <- "results/07_immune_subclustering/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ── 1. Load annotated full Seurat object ─────────────────────────────────────
seurat_merged <- readRDS(file.path(input_dir, "seurat_merged_clustered.rds"))

## ── 2. Subset Immune Cells ───────────────────────────────────────────────────
immune <- subset(seurat_merged, idents = "Immune Cells")
DefaultAssay(immune) <- "RNA"
immune[["RNA"]] <- JoinLayers(immune[["RNA"]])

cat("Immune cells total:", ncol(immune), "\n")
print(table(immune$treatment))

## ── 3. Independent re-clustering of immune cells ─────────────────────────────
# Re-normalize and scale within the immune subset to capture fine-grained
# immune heterogeneity lost in the full-atlas clustering.
immune <- NormalizeData(immune)
immune <- FindVariableFeatures(immune, selection.method = "vst", nfeatures = 2000)
immune <- ScaleData(immune)
immune <- RunPCA(immune, npcs = 30)

# Elbow plot saved for documentation
p_elbow <- ElbowPlot(immune, ndims = 30)
ggsave(file.path(output_dir, "immune_elbow_plot.pdf"), p_elbow, width = 5, height = 4)

# 15 PCs, resolution sweep – res = 0.3 selected (8 biologically coherent clusters)
immune <- FindNeighbors(immune, dims = 1:15)
immune <- FindClusters(immune, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))
immune <- RunUMAP(immune, dims = 1:15)

cat("\nCluster sizes at res = 0.3:\n")
print(table(immune$originalexp_snn_res.0.3))

# Set working identity to the selected resolution
immune$seurat_clusters <- immune$originalexp_snn_res.0.3
Idents(immune) <- "seurat_clusters"

## ── 4. Marker genes for subcluster annotation ─────────────────────────────────
markers_immune <- FindAllMarkers(
    immune,
    assay           = "RNA",
    only.pos        = TRUE,
    test.use        = "wilcox",
    min.pct         = 0.25,
    logfc.threshold = 0.25
)

top10_immune <- markers_immune %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 10)

write.csv(markers_immune,
    file.path(output_dir, "markers_immune_subclusters.csv"), row.names = FALSE)

## ── 5. Cell-type annotation ───────────────────────────────────────────────────
# Cluster labels assigned based on canonical immune marker expression.
# Cluster 0 (Neurons contamination) removed after annotation.
immune <- RenameIdents(
    immune,
    "0" = "Neurons",           # contaminating neurons – will be removed
    "1" = "ITRM",              # Immunosuppressive tissue-resident macrophages
    "2" = "Tregs",             # Regulatory T cells
    "3" = "CD8+ T",            # CD8+ T cells
    "4" = "DC",                # Dendritic cells
    "5" = "B cells",           # B cells
    "6" = "Act. DC",           # Activated dendritic cells
    "7" = "Diff. neurons",     # Differentiated neurons (contamination)
    "8" = "pDC",               # Plasmacytoid dendritic cells
    "9" = "Mast cells"         # Mast cells
)

immune$immune_annotation <- factor(
    as.character(Idents(immune)),
    levels = sort(unique(as.character(Idents(immune))))
)

cat("\nAnnotation × treatment:\n")
print(table(immune$immune_annotation, immune$treatment))

# Remove neuronal contaminants (clusters 0)
immune <- subset(immune, subset = immune_annotation %in%
    c("Neurons") == FALSE)
immune$immune_annotation <- droplevels(immune$immune_annotation)

cat("\nFinal immune subtypes:\n")
print(table(immune$immune_annotation, immune$treatment))

saveRDS(immune, file.path(output_dir, "NK_immune.rds"))

## ── 6. ON vs OFF differential expression ─────────────────────────────────────
Idents(immune) <- "treatment"

# All immune cells combined
deg_all <- FindMarkers(immune, ident.1 = "ON", ident.2 = "OFF",
    assay = "RNA", test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(deg_all, file.path(output_dir, "deg_ON_OFF_immune_all.csv"),
    row.names = TRUE)

# Subcluster 1 (ITRM) specifically
immune_sub1 <- subset(immune,
    subset = originalexp_snn_res.0.3 == "1")
Idents(immune_sub1) <- "treatment"
deg_sub1 <- FindMarkers(immune_sub1, ident.1 = "ON", ident.2 = "OFF",
    assay = "RNA", test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(deg_sub1, file.path(output_dir, "deg_ON_OFF_immune_subcluster1.csv"),
    row.names = TRUE)

# Fisher's exact test: cluster composition ON vs OFF
tab <- table(immune$originalexp_snn_res.0.3, immune$treatment)
fisher_res <- fisher.test(tab, simulate.p.value = TRUE)
cat(sprintf("\nFisher's exact test for immune cluster × treatment: p = %.4f\n",
    fisher_res$p.value))

## ── 7. Color palette and short labels ────────────────────────────────────────
cell_colors <- c(
    "Act. DC"      = "#E69F00",
    "B cells"      = "#56B4E9",
    "CD8+ T"       = "#882255",
    "DC"           = "#009E73",
    "Diff. neurons"= "#0072B2",
    "ITRM"         = "#D55E00",
    "Mast cells"   = "#CC79A7",
    "pDC"          = "#F0E442",
    "Tregs"        = "#117733"
)

immune$label_short <- as.character(immune$immune_annotation)

## ── 8. UMAP panels (Fig 6D) ───────────────────────────────────────────────────
# Fig 6D: Split UMAP – ON panel (labeled centroids), OFF panel (larger points,
#          sparse due to dramatic immune contraction upon H3.3K27M withdrawal)

get_centroids <- function(seurat_obj, subset_cells = NULL) {
    emb <- as.data.frame(Embeddings(seurat_obj, "umap"))
    colnames(emb) <- c("UMAP_1", "UMAP_2")
    emb$label <- seurat_obj$label_short
    if (!is.null(subset_cells)) emb <- emb[subset_cells, ]
    emb %>%
        group_by(label) %>%
        summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = "drop")
}

all_emb <- as.data.frame(Embeddings(immune, "umap"))
colnames(all_emb) <- c("UMAP_1", "UMAP_2")
all_emb$annotation <- immune$label_short

on_cells  <- WhichCells(immune, expression = treatment == "ON")
off_cells <- WhichCells(immune, expression = treatment == "OFF")
on_emb    <- all_emb[on_cells, ]
off_emb   <- all_emb[off_cells, ]

umap_theme <- theme_classic(base_size = 11) +
    theme(
        legend.position  = "none",
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        axis.title       = element_text(size = 9, color = "grey40"),
        plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

# ON panel with repelled centroid labels
centroids_on <- get_centroids(immune, on_cells)
p_on <- ggplot() +
    geom_point(data = all_emb, aes(UMAP_1, UMAP_2),
        color = "grey88", size = 0.6, alpha = 0.5) +
    geom_point(data = on_emb, aes(UMAP_1, UMAP_2, color = annotation),
        size = 1.2, alpha = 0.8) +
    geom_label_repel(data = centroids_on,
        aes(UMAP_1, UMAP_2, label = label, color = label),
        size = 2.8, fontface = "bold",
        fill = alpha("white", 0.75),
        label.size = 0.25, label.padding = unit(0.15, "lines"),
        max.overlaps = Inf, force = 3, seed = 42, show.legend = FALSE) +
    scale_color_manual(values = cell_colors) +
    labs(title = "ON", x = "UMAP 1", y = "UMAP 2") +
    umap_theme

# OFF panel – sparse cells shown as larger points, no labels
p_off <- ggplot() +
    geom_point(data = all_emb, aes(UMAP_1, UMAP_2),
        color = "grey88", size = 0.6, alpha = 0.5) +
    geom_point(data = off_emb,
        aes(UMAP_1, UMAP_2, color = annotation, fill = annotation),
        size = 2.5, alpha = 0.9, stroke = 0.3, shape = 21) +
    scale_color_manual(values = cell_colors) +
    scale_fill_manual(values  = cell_colors) +
    labs(title = "OFF", x = "UMAP 1", y = "") +
    umap_theme

# Shared legend
legend_plot <- ggplot(all_emb, aes(UMAP_1, UMAP_2, color = annotation)) +
    geom_point(size = 2) +
    scale_color_manual(values = cell_colors, name = NULL) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
    theme(legend.text = element_text(size = 8),
          legend.key.size = unit(0.5, "cm"))
shared_legend <- cowplot::get_legend(legend_plot)

fig6D <- cowplot::plot_grid(
    p_on, p_off, shared_legend,
    nrow = 1, rel_widths = c(1, 1, 0.55)
)

ggsave(file.path(output_dir, "immune_UMAP_ON_OFF_Fig6D.png"),
    fig6D, width = 9, height = 4, dpi = 300)

## ── 9. Composition barplot (Fig 6E) ──────────────────────────────────────────
# Fig 6E: Stacked barplot of immune subtype proportions, ON vs OFF
bar_df <- data.frame(
    annotation = immune$label_short,
    treatment  = immune$treatment
) %>%
    count(treatment, annotation) %>%
    group_by(treatment) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup() %>%
    mutate(
        treatment  = factor(treatment, levels = c("ON", "OFF")),
        annotation = factor(annotation, levels = rev(names(cell_colors)))
    )

fig6E <- ggplot(bar_df, aes(x = treatment, y = pct, fill = annotation)) +
    geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = cell_colors, name = NULL) +
    scale_y_continuous(
        expand = c(0, 0), limits = c(0, 101),
        breaks = seq(0, 100, 25),
        labels = label_percent(scale = 1)
    ) +
    labs(x = NULL, y = "Cell type composition (%)") +
    theme_classic(base_size = 11) +
    theme(
        axis.text.x      = element_text(size = 11, color = "black", face = "bold"),
        axis.text.y      = element_text(size = 9,  color = "black"),
        axis.title.y     = element_text(size = 10, color = "black"),
        axis.ticks.x     = element_blank(),
        axis.line        = element_line(color = "black", linewidth = 0.4),
        legend.text      = element_text(size = 8),
        legend.key.size  = unit(0.4, "cm"),
        legend.position  = "right",
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        plot.margin      = margin(10, 10, 10, 10)
    ) +
    guides(fill = guide_legend(reverse = TRUE, ncol = 1))

ggsave(file.path(output_dir, "immune_composition_barplot_Fig6E.png"),
    fig6E, width = 5, height = 4.5, dpi = 300)

cat("\nImmune subclustering complete. Outputs saved to:", output_dir, "\n")
