##=============================================================================
## 03_differential_expression.R
##
## Project : iH3.3K27M DMG scRNA-seq (Khairkhah et al.)
## Purpose : (A) Cluster marker genes (FindAllMarkers) used for cell-type
##           annotation; (B) ON vs OFF differential expression within
##           individual cell-type clusters; (C) Stemness/differentiation
##           heatmaps (Sup Fig 7F); (D) H3f3a expression violin plot
##           comparing ON vs OFF (Sup Fig 7B); (E) Glutamate/aspartate
##           transporter and receptor violin plots in DMG tumor cells
##           (Slc1a3, Glud1, Got1, Got2, Gria3, Gria4) (Sup Fig 7E).
## Methods : Seurat v5.1.0 Wilcoxon rank-sum test on the non-integrated RNA
##           assay. Bonferroni-adjusted p-values.
## Utilities: utilities/violin_plot_celltype.R
##            – plot_violin_panels(): wraps VlnPlot with adj.P and log2FC
##              annotations per gene panel (Sup Fig 7D/E/G)
##            – make_violin_plot_list(): returns a list of ggplot objects
##              for finer layout control
## Input   : seurat_merged_clustered.rds
## Output  : DE tables (xlsx) + figures
##=============================================================================

## ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(openxlsx)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(patchwork)    # required by violin_plot_celltype.R

## ── Utility functions ─────────────────────────────────────────────────────────
# violin_plot_celltype.R provides:
#   plot_violin_panels()      – combined patchwork figure of per-gene violins
#   make_violin_plot_list()   – list of individual ggplot objects
#   make_de_label() / fmt_adj_p() – internal label helpers
source("utilities/violin_plot_celltype.R")

## ── Paths ─────────────────────────────────────────────────────────────────────
input_dir  <- "results/02_clustering/"
output_dir <- "results/03_DE/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ── Load ──────────────────────────────────────────────────────────────────────
seurat_merged <- readRDS(file.path(input_dir, "seurat_merged_clustered.rds"))

# Ensure DE is run on the original (non-integrated) RNA assay
DefaultAssay(seurat_merged) <- "RNA"
seurat_merged <- JoinLayers(seurat_merged)   # collapse split layers for FindMarkers

## ══════════════════════════════════════════════════════════════════════════════
## A. Cluster marker genes (used for cell-type annotation)
## ══════════════════════════════════════════════════════════════════════════════
Idents(seurat_merged) <- "CellType_alpha"

de_clusters <- FindAllMarkers(
    seurat_merged,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.1,
    test.use        = "wilcox"
)
de_clusters_filt <- de_clusters %>% filter(p_val_adj < 0.05)

# Top 5 markers per cluster (for dot plot)
de_top5 <- de_clusters_filt %>%
    group_by(cluster) %>%
    arrange(p_val_adj, desc(avg_log2FC)) %>%
    slice_head(n = 5) %>%
    ungroup()

# Top 50 for supplementary table
de_top50 <- de_clusters_filt %>%
    group_by(cluster) %>%
    arrange(p_val_adj) %>%
    slice_head(n = 50) %>%
    ungroup()

write.xlsx(
    list("All_markers"  = de_clusters_filt,
         "Top50_per_cluster" = de_top50),
    file.path(output_dir, "cluster_markers.xlsx")
)

# Dot plot – top 5 markers per cluster
p_dot <- DotPlot(seurat_merged, features = unique(de_top5$gene)) +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(title = "Top 5 marker genes per cell-type cluster")

ggsave(file.path(output_dir, "dotplot_top5_markers.pdf"),
       p_dot, width = 14, height = 18)

## ══════════════════════════════════════════════════════════════════════════════
## B. ON vs OFF DE – per cluster
##    (Sup Fig 7D/E: excitatory neuron cluster 3, inhibitory neuron cluster 20)
## ══════════════════════════════════════════════════════════════════════════════
clusters_of_interest <- c("Excitatory Neurons", "Inhibitory Neurons")

de_list <- lapply(clusters_of_interest, function(ct) {
    sub_obj <- subset(seurat_merged, subset = CellType_alpha == ct)
    Idents(sub_obj) <- "treatment"
    FindMarkers(sub_obj, ident.1 = "ON", ident.2 = "OFF",
                test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0)
})
names(de_list) <- clusters_of_interest

write.xlsx(de_list, file.path(output_dir, "DE_ON_vs_OFF_per_celltype.xlsx"),
           rowNames = TRUE)

# ── Violin plots using utilities/violin_plot_celltype.R ──────────────────────
# plot_violin_panels() handles:
#   • VlnPlot generation per gene
#   • adj.P and log2FC annotation from the DE table
#   • patchwork layout into a single figure
#
# The de_table argument must be a FindMarkers() result with genes as rownames.

# Excitatory neuron markers (Sup Fig 7D)
excit_markers <- c("Gabra6", "Cbln1", "Cbln3", "Grm4", "Grin2c", "Neurod1",
                   "Syt2", "Adcy1")
excit_sub <- subset(seurat_merged, subset = CellType_alpha == "Excitatory Neurons")
excit_sub$treatment <- factor(excit_sub$treatment, levels = c("ON", "OFF"))

p_excit <- plot_violin_panels(
    seu       = excit_sub,
    genes     = excit_markers,
    de_table  = de_list[["Excitatory Neurons"]],
    group.by  = "treatment",
    cols      = c("ON" = "#F1C40F", "OFF" = "#A7A9AC"),
    y_lim     = c(0, 7),
    ncol      = 4
)
ggsave(file.path(output_dir, "violin_ExcitatoryNeuron_markers_SupFig7D.pdf"),
       p_excit, width = 14, height = 6)

# Inhibitory neuron markers (Sup Fig 7E)
inhib_markers <- c("Cnr1", "Tfap2b", "Pvalb")
inhib_sub <- subset(seurat_merged, subset = CellType_alpha == "Inhibitory Neurons")
inhib_sub$treatment <- factor(inhib_sub$treatment, levels = c("ON", "OFF"))

p_inhib <- plot_violin_panels(
    seu       = inhib_sub,
    genes     = inhib_markers,
    de_table  = de_list[["Inhibitory Neurons"]],
    group.by  = "treatment",
    cols      = c("ON" = "#F1C40F", "OFF" = "#A7A9AC"),
    y_lim     = c(0, 7),
    ncol      = 3
)
ggsave(file.path(output_dir, "violin_InhibitoryNeuron_markers_SupFig7E.pdf"),
       p_inhib, width = 10, height = 4)

## ══════════════════════════════════════════════════════════════════════════════
## C. Stemness / Differentiation heatmap (Sup Fig 7F)
##    Shows ON vs OFF average log2FC for curated gene sets in DMG tumor cells
## ══════════════════════════════════════════════════════════════════════════════
stemness_genes <- c("Sox2", "Pax6", "Ascl1", "Hes1", "Hes5", "Prom1",
                    "Itga6", "Lgr5", "Bmi1", "Pcna")
diff_genes     <- c("Aqp4", "Rbfox3", "Pdgfra", "Sox10", "Cnp", "Mbp",
                    "Plp1", "Mog", "Mag")

tumor_cells <- subset(seurat_merged, subset = CellType_alpha == "DMG Tumor")
Idents(tumor_cells) <- "treatment"

all_genes <- c(stemness_genes, diff_genes)
de_tumor <- FindMarkers(tumor_cells, ident.1 = "ON", ident.2 = "OFF",
                        features = all_genes, logfc.threshold = 0,
                        min.pct = 0, test.use = "wilcox")

# Z-score the log2FC values across genes for heatmap
mat <- matrix(de_tumor$avg_log2FC, ncol = 1,
              dimnames = list(rownames(de_tumor), "ON_vs_OFF_log2FC"))
mat_z <- scale(mat)

row_ann <- rowAnnotation(
    GeneSet = ifelse(rownames(mat_z) %in% stemness_genes,
                     "Stemness", "Differentiation"),
    col = list(GeneSet = c("Stemness" = "#E41A1C", "Differentiation" = "#377EB8")),
    show_annotation_name = FALSE
)

ht <- Heatmap(
    mat_z,
    name              = "Z score",
    col               = colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C")),
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    left_annotation   = row_ann,
    column_title      = "ON vs OFF (DMG Tumor)",
    row_names_gp      = gpar(fontsize = 9),
    column_names_gp   = gpar(fontsize = 9)
)

pdf(file.path(output_dir, "heatmap_stemness_diff_SupFig7F.pdf"),
    width = 3.5, height = 6)
draw(ht)
dev.off()

## ══════════════════════════════════════════════════════════════════════════════
## D. H3f3a expression violin – ON vs OFF (Sup Fig 7B)
## ══════════════════════════════════════════════════════════════════════════════
df_h3 <- FetchData(seurat_merged, vars = c("H3f3a", "treatment"))
df_h3$treatment <- factor(df_h3$treatment, levels = c("ON", "OFF"))

mypal <- pal_npg("nrc", alpha = 1)(2)

p_h3 <- ggviolin(df_h3, x = "treatment", y = "H3f3a", fill = "treatment",
    palette  = mypal,
    add      = "boxplot",
    add.params = list(fill = "white")) +
    stat_compare_means(comparisons = list(c("ON", "OFF")),
                       method = "wilcox.test") +
    labs(x = "Treatment", y = "H3f3a normalized expression") +
    NoLegend()

ggsave(file.path(output_dir, "violin_H3f3a_ON_OFF_SupFig7B.pdf"),
       p_h3, width = 3.5, height = 4.5)

# Wilcoxon test summary
wt <- wilcox.test(H3f3a ~ treatment, data = df_h3, exact = FALSE)
cat(sprintf("\nH3f3a ON vs OFF Wilcoxon p-value: %.3e\n", wt$p.value))

## ══════════════════════════════════════════════════════════════════════════════
## E. Glutamate/aspartate transporter & receptor violins – DMG tumor cells
##    (Sup Fig 7E)
##
##    Genes: Slc1a3 (GLAST – glial glutamate transporter), Glud1 (glutamate
##    dehydrogenase), Got1/Got2 (glutamate-oxaloacetate transaminases 1/2),
##    Gria3 (GluA3 – AMPAR subunit), Gria4 (GluA4 – AMPAR subunit).
##    These markers capture glutamate uptake and AMPAR-mediated signaling
##    within DMG tumor cells (clusters 16+17) ON vs OFF.
## ══════════════════════════════════════════════════════════════════════════════

# Subset DMG tumor cells (clusters 16 and 17 = PDGFRA+ OPC Early/Differentiating)
dmg_tumor <- subset(seurat_merged, subset = CellType_alpha == "DMG Tumor")
dmg_tumor$treatment <- factor(dmg_tumor$treatment, levels = c("ON", "OFF"))

# ON vs OFF DE within DMG tumor cells (all genes, for annotation)
Idents(dmg_tumor) <- "treatment"
de_dmg <- FindMarkers(dmg_tumor, ident.1 = "ON", ident.2 = "OFF",
    test.use = "wilcox", min.pct = 0, logfc.threshold = 0)

write.xlsx(de_dmg,
    file.path(output_dir, "DE_ON_vs_OFF_DMG_tumor.xlsx"), rowNames = TRUE)

# Violin panels using utilities/violin_plot_celltype.R
# plot_violin_panels() pulls adj.P and log2FC directly from the DE table
# and annotates each gene's panel automatically.
glut_genes <- c("Slc1a3", "Glud1", "Got1", "Got2", "Gria3", "Gria4")

p_dmg_glut <- plot_violin_panels(
    seu      = dmg_tumor,
    genes    = glut_genes,
    de_table = de_dmg,
    group.by = "treatment",
    cols     = c("ON" = "#F1C40F", "OFF" = "#A7A9AC"),
    y_lim    = c(0, 7),
    ncol     = 3
)

ggsave(file.path(output_dir, "violin_DMG_glutamate_markers_SupFig7E.pdf"),
    p_dmg_glut, width = 9, height = 6)
