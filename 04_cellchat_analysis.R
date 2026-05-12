##=============================================================================
## 04_cellchat_analysis.R
##
## Project : iH3.3K27M DMG scRNA-seq (Khairkhah et al.)
## Purpose : Build CellChat objects for ON and OFF conditions separately,
##           compute communication probabilities, centrality scores, and
##           communication patterns; merge objects for comparative analysis.
## Methods : CellChat v1.6.1, CellChatDB.mouse ligand-receptor database,
##           truncated mean method.
## Utilities: utilities/rankNet_patched.R
##            – rankNet_patched(): replacement for CellChat::rankNet() that
##              (a) returns the ggplot data frame when return.data = TRUE,
##              (b) colors ON/OFF bars as gold/grey, and
##              (c) removes non-significant pathways cleanly.
##            utilities/netAnalysis_signalingRole_heatmap2.R
##            – netAnalysis_signalingRole_heatmap2(): extended heatmap that
##              accepts named color.use vectors and fixed axis limits
##              (top.bar.ylim, right.bar.ylim) for cross-condition comparability.
## Input   : seurat_merged_clustered.rds
## Output  : cellchaton_analysed.rds, cellchatoff_analysed.rds,
##           cellchat_on_off_merged.rds
## Figures : Fig 5B (rankNet bar), Fig 5D (overall signaling heatmap),
##           Fig 6A (circle – number of interactions),
##           Fig 6B (circle – interaction strength),
##           Fig 6G (outgoing/incoming signaling heatmaps),
##           Sup Fig 7H (differential circles), Sup Fig 7I (bar totals)
##=============================================================================

## ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(CellChat)    # v1.6.1
library(ggplot2)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)

## ── Utility functions ─────────────────────────────────────────────────────────
# rankNet_patched.R: drop-in replacement for CellChat::rankNet() in
# comparison mode. Key additions vs the original:
#   • return.data = TRUE returns list(signaling.contribution = df, gg.obj = gg)
#   • ON/OFF bars colored gold (#F1C40F) and grey (#A7A9AC) by default
#   • Non-significant pathways suppressed rather than shown in black
source("utilities/rankNet_patched.R")

# netAnalysis_signalingRole_heatmap2.R: extended version of CellChat's
# netAnalysis_signalingRole_heatmap() that accepts:
#   • Named color.use vector keyed by cell-type name
#   • top.bar.ylim and right.bar.ylim for fixed axis scales across conditions
source("utilities/netAnalysis_signalingRole_heatmap2.R")

## ── Paths ─────────────────────────────────────────────────────────────────────
input_dir  <- "results/02_clustering/"
output_dir <- "results/04_cellchat/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ── Helper: build one CellChat object from a Seurat subset ───────────────────
build_cellchat <- function(seurat_obj, condition, db = CellChatDB.mouse) {
    sub  <- subset(seurat_obj, subset = treatment == condition)
    data <- GetAssayData(sub, assay = "RNA", layer = "data")
    meta <- sub@meta.data[, "CellType_alpha", drop = FALSE]
    
    cc <- createCellChat(object = data, meta = meta, group.by = "CellType_alpha")
    cc@DB <- db
    
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1)
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
    
    cat(sprintf("[%s] Pathways inferred: %d\n", condition,
                length(cc@netP$pathways)))
    cc
}

## ── 1. Load Seurat object ─────────────────────────────────────────────────────
seurat_merged <- readRDS(file.path(input_dir, "seurat_merged_clustered.rds"))
DefaultAssay(seurat_merged) <- "RNA"
seurat_merged <- JoinLayers(seurat_merged)

## ── 2. Build separate CellChat objects ───────────────────────────────────────
# Use CellChatDB.mouse – curated mouse ligand-receptor interactions
CellChatDB <- CellChatDB.mouse

cellchaton  <- build_cellchat(seurat_merged, "ON",  db = CellChatDB)
cellchatoff <- build_cellchat(seurat_merged, "OFF", db = CellChatDB)

saveRDS(cellchaton,  file.path(output_dir, "cellchaton_analysed.rds"))
saveRDS(cellchatoff, file.path(output_dir, "cellchatoff_analysed.rds"))

## ── 3. Communication patterns (selectK, identifyCommunicationPatterns) ────────
selectK(cellchaton, pattern = "outgoing")
cellchaton <- identifyCommunicationPatterns(cellchaton, pattern = "outgoing",
    k = 2, width = 5, height = 15)
selectK(cellchaton, pattern = "incoming")
cellchaton <- identifyCommunicationPatterns(cellchaton, pattern = "incoming",
    k = 2, width = 5, height = 15)

selectK(cellchatoff, pattern = "outgoing")
cellchatoff <- identifyCommunicationPatterns(cellchatoff, pattern = "outgoing",
    k = 2, width = 5, height = 15)
selectK(cellchatoff, pattern = "incoming")
cellchatoff <- identifyCommunicationPatterns(cellchatoff, pattern = "incoming",
    k = 4, width = 5, height = 15)

## ── 4. Merge for comparative analysis ────────────────────────────────────────
object.list <- list(ON = cellchaton, OFF = cellchatoff)
cellchat    <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat, file.path(output_dir, "cellchat_on_off_merged.rds"))

## ── 5. Pathway information flow bar chart (Fig 5B) ───────────────────────────
# Fig 5B: Stacked bar chart of relative information flow across all pathways

res <- rankNet_patched(cellchat, mode = "comparison", stacked = TRUE,
                       do.stat = TRUE, return.data = TRUE)

p_rank <- res$gg.obj +
    xlab("Cell-chat signaling pathways") +
    theme_classic(base_size = 12) +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))

ggsave(file.path(output_dir, "rankNet_stacked_Fig5B.pdf"),
       p_rank, width = 4, height = 6)

# Export the underlying data for transparency
write.csv(res$signaling.contribution,
          file.path(output_dir, "rankNet_data_Fig5B.csv"), row.names = FALSE)

## ── 6. Overall signaling pattern heatmaps (Fig 5D) ───────────────────────────
# Fig 5D: Heatmap of overall signaling role (sender/receiver) per cell type
# for a curated subset of pathways shown in ON and OFF panels.

mycol <- setNames(
    c("#AEC6CF","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
      "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#A6CEE3"),
    c("Astrocyte","DMG Tumor","Endothelial Cells","Ependymal Cells",
      "Epithelial Cells","Excitatory Neurons","Fibroblasts","Immune Cells",
      "Inhibitory Neurons","Microglia","Neurons","Oligodendrocytes",
      "Purkinje Neurons")
)

all_pathways    <- union(object.list[[1]]@netP$pathways,
                          object.list[[2]]@netP$pathways)
subset_pathways <- c("NRXN","NCAM","CADM","PTN","CNTN","PSAP","LAMININ","APP",
                     "EPHA","NT","OPIOID","SOMATOSTATIN","NPY","EDN",
                     "CD200","CD45","CEACAM")

ht1 <- netAnalysis_signalingRole_heatmap2(object.list[[1]], pattern = "all",
    signaling = subset_pathways, title = "ON",
    width = 5, height = 8, color.heatmap = "OrRd",
    top.bar.ylim = c(0, 1), right.bar.ylim = c(0, 0.8), color.use = mycol)
ht2 <- netAnalysis_signalingRole_heatmap2(object.list[[2]], pattern = "all",
    signaling = subset_pathways, title = "OFF",
    width = 5, height = 8, color.heatmap = "OrRd",
    top.bar.ylim = c(0, 1), right.bar.ylim = c(0, 0.8), color.use = mycol)

pdf(file.path(output_dir, "overall_signaling_heatmap_Fig5D.pdf"),
    width = 10, height = 7)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

# Outgoing-only heatmap (Fig 6G top)
ht1_out <- netAnalysis_signalingRole_heatmap2(object.list[[1]], pattern = "outgoing",
    signaling = subset_pathways, title = "ON",
    width = 5, height = 8, color.heatmap = "OrRd",
    top.bar.ylim = c(0,1), right.bar.ylim = c(0, 0.8), color.use = mycol)
ht2_out <- netAnalysis_signalingRole_heatmap2(object.list[[2]], pattern = "outgoing",
    signaling = subset_pathways, title = "OFF",
    width = 5, height = 8, color.heatmap = "OrRd",
    top.bar.ylim = c(0,1), right.bar.ylim = c(0, 0.8), color.use = mycol)

pdf(file.path(output_dir, "outgoing_signaling_heatmap_Fig6G_top.pdf"),
    width = 10, height = 7)
draw(ht1_out + ht2_out, ht_gap = unit(0.5, "cm"))
dev.off()

# Incoming-only heatmap (Fig 6G bottom)
ht1_in <- netAnalysis_signalingRole_heatmap2(object.list[[1]], pattern = "incoming",
    signaling = subset_pathways, title = "ON",
    width = 5, height = 8, color.heatmap = "OrRd",
    top.bar.ylim = c(0,1), right.bar.ylim = c(0, 0.8), color.use = mycol)
ht2_in <- netAnalysis_signalingRole_heatmap2(object.list[[2]], pattern = "incoming",
    signaling = subset_pathways, title = "OFF",
    width = 5, height = 8, color.heatmap = "OrRd",
    top.bar.ylim = c(0,1), right.bar.ylim = c(0, 0.8), color.use = mycol)

pdf(file.path(output_dir, "incoming_signaling_heatmap_Fig6G_bottom.pdf"),
    width = 10, height = 7)
draw(ht1_in + ht2_in, ht_gap = unit(0.5, "cm"))
dev.off()

## ── 7. Differential interactions: number & strength (Fig 6A, 6B / Sup 7H) ───
# Fig 6A: Number of interactions circle plot (ON vs OFF)
# Fig 6B: Interaction strength circle plot (ON vs OFF)

weight_max <- max(
    sapply(object.list, function(x) max(x@net$weight, na.rm = TRUE))
)

png(file.path(output_dir, "circle_nInteractions_Fig6A.png"),
    width = 12, height = 6, units = "in", res = 600, bg = "white")
par(mfrow = c(1, 2), xpd = TRUE)
for (i in seq_along(object.list)) {
    netVisual_circle(object.list[[i]]@net$count,
        vertex.weight     = as.numeric(table(object.list[[i]]@idents)),
        weight.scale      = TRUE,
        label.edge        = FALSE,
        title.name        = paste("Number of interactions –", names(object.list)[i]),
        color.use         = mycol)
}
dev.off()

png(file.path(output_dir, "circle_interactionStrength_Fig6B.png"),
    width = 12, height = 6, units = "in", res = 600, bg = "white")
par(mfrow = c(1, 2), xpd = TRUE)
for (i in seq_along(object.list)) {
    netVisual_circle(object.list[[i]]@net$weight,
        vertex.weight     = as.numeric(table(object.list[[i]]@idents)),
        weight.scale      = TRUE,
        label.edge        = FALSE,
        title.name        = paste("Interaction strength –", names(object.list)[i]),
        color.use         = mycol)
}
dev.off()

# Differential circle plots (Sup Fig 7H)
png(file.path(output_dir, "circle_diff_interactions_SupFig7H.png"),
    width = 12, height = 6, units = "in", res = 600, bg = "white")
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat, comparison = c(1, 2), measure = "count",
    weight.scale = TRUE, arrow.size = 0.1)
netVisual_diffInteraction(cellchat, comparison = c(1, 2), measure = "weight",
    weight.scale = TRUE, arrow.size = 0.1)
dev.off()

cat("\nCellChat analysis complete. Outputs saved to:", output_dir, "\n")
