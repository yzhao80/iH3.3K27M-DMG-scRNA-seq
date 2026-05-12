##=============================================================================
## 06_cellchat_heatmap_pathways.R
##
## Project : iH3.3K27M DMG scRNA-seq (Khairkhah et al.)
## Purpose : Pathway-level heatmaps comparing signaling strength and
##           ligand-receptor interactions between ON and OFF for key pathways.
##           Also generates the diverging bar chart of opposing signaling
##           (Fig 5C) and the LR-pair dot/heatmap plots (Fig 5G, Fig 6F).
## Utilities: utilities/rankNet_patched.R
##            – rankNet_patched(): used to extract signaling contribution data
##              (return.data = TRUE) for the diverging bar construction.
##            utilities/netAnalysis_signalingRole_heatmap2.R
##            – netAnalysis_signalingRole_heatmap2(): fixed-scale heatmap for
##              outgoing/incoming signaling patterns across ON and OFF.
## Figures : Fig 5C (diverging bar – DMG tumor cell pathways)
##           Fig 5G (NCAM/NRXN/PTN/CADM LR heatmap)
##           Fig 5H (LR dot plot – ON vs OFF)
##           Fig 6F (EPHA/NT/OPIOID/SOMATOSTATIN LR dot plot)
## Input   : cellchaton_analysed.rds, cellchatoff_analysed.rds,
##           cellchat_on_off_merged.rds, Net_on.txt, Net_off.txt
##=============================================================================

## ── Libraries ─────────────────────────────────────────────────────────────────
library(CellChat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

## ── Utility functions ─────────────────────────────────────────────────────────
source("utilities/rankNet_patched.R")
source("utilities/netAnalysis_signalingRole_heatmap2.R")

## ── Paths ─────────────────────────────────────────────────────────────────────
input_dir  <- "results/04_cellchat/"
output_dir <- "results/06_heatmap_pathways/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ── Load CellChat objects ─────────────────────────────────────────────────────
cellchaton  <- readRDS(file.path(input_dir, "cellchaton_analysed.rds"))
cellchatoff <- readRDS(file.path(input_dir, "cellchatoff_analysed.rds"))
cellchat    <- readRDS(file.path(input_dir, "cellchat_on_off_merged.rds"))
object.list <- list(ON = cellchaton, OFF = cellchatoff)

## ── Internal helper: extract signaling changes for one cell type ──────────────
# Mirrors CellChat::netAnalysis_signalingChanges_scatter logic but returns a
# tidy data frame for downstream use.
extract_signalingChanges_df <- function(cc_merged, idents.use,
                                        comparison = c(1, 2)) {
    dataset.name <- names(cc_merged@object.list)
    
    mat.all <- lapply(cc_merged@object.list, function(obj) {
        prob  <- obj@netP$prob
        out   <- apply(prob[idents.use, , , drop = FALSE], 3,
                       function(x) sum(x, na.rm = TRUE))
        inc   <- apply(prob[, idents.use, , drop = FALSE], 3,
                       function(x) sum(x, na.rm = TRUE))
        cbind(outgoing = out, incoming = inc)
    })
    
    all_pw   <- union(rownames(mat.all[[1]]), rownames(mat.all[[2]]))
    blank    <- matrix(0, nrow = length(all_pw), ncol = 2,
                       dimnames = list(all_pw, c("outgoing", "incoming")))
    mat.all.merged.use <- lapply(mat.all, function(m) {
        r          <- blank
        r[rownames(m), ] <- m
        r
    })
    
    mat.diff <- mat.all.merged.use[[comparison[2]]] -
                mat.all.merged.use[[comparison[1]]]
    keep     <- rowSums(abs(mat.diff)) > 0
    mat.diff <- mat.diff[keep, , drop = FALSE]
    
    idx.specific <- mat.all.merged.use[[1]] * mat.all.merged.use[[2]]
    mat.sum      <- mat.all.merged.use[[2]] + mat.all.merged.use[[1]]
    
    out.sp <- rownames(idx.specific)[
        rownames(idx.specific) %in% rownames(mat.diff) &
        (mat.sum[rownames(mat.diff), 1] != 0) &
        (idx.specific[rownames(mat.diff), 1] == 0)]
    in.sp  <- rownames(idx.specific)[
        rownames(idx.specific) %in% rownames(mat.diff) &
        (mat.sum[rownames(mat.diff), 2] != 0) &
        (idx.specific[rownames(mat.diff), 2] == 0)]
    
    df            <- as.data.frame(mat.diff)
    df$specificity <- "Shared"
    df$specificity[rownames(df) %in% out.sp & rowSums(mat.diff >= 0) == 2] <-
        paste0(dataset.name[comparison[2]], " specific")
    df$specificity[rownames(df) %in% in.sp & rowSums(mat.diff <= 0) == 2] <-
        paste0(dataset.name[comparison[1]], " specific")
    df$labels <- rownames(df)
    df
}

## ══════════════════════════════════════════════════════════════════════════════
## A. Diverging bar chart – opposing signaling in DMG Tumor Cells (Fig 5C)
## ══════════════════════════════════════════════════════════════════════════════
get_pathway_strength <- function(cc_obj) {
    probs <- cc_obj@netP$prob
    data.frame(
        pathway  = dimnames(probs)[[3]],
        strength = apply(probs, 3, function(x) sum(x, na.rm = TRUE)),
        stringsAsFactors = FALSE
    )
}

df_on  <- get_pathway_strength(cellchaton)  %>% rename(ON  = strength)
df_off <- get_pathway_strength(cellchatoff) %>% rename(OFF = strength)

df_strength <- full_join(df_on, df_off, by = "pathway") %>%
    mutate(ON = replace_na(ON, 0), OFF = replace_na(OFF, 0))

df.scatter      <- extract_signalingChanges_df(cellchat, "DMG Tumor")
off_specific_pw <- df.scatter %>%
    filter(specificity == "OFF specific") %>%
    pull(labels)

must_keep   <- c("NCAM", "PTN", "NRXN", "CADM")
pathways_use <- union(must_keep, off_specific_pw)

plot_df <- df_strength %>%
    filter(pathway %in% pathways_use) %>%
    mutate(max_strength = pmax(ON, OFF)) %>%
    arrange(desc(max_strength)) %>%
    mutate(pathway = factor(pathway, levels = rev(pathway)))

plot_long <- plot_df %>%
    pivot_longer(c("ON", "OFF"), names_to = "condition", values_to = "strength") %>%
    mutate(plot_value = ifelse(condition == "ON", -strength, strength))

p_div <- ggplot(plot_long,
                aes(x = pathway, y = plot_value, fill = condition)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c("ON" = "#F1C40F", "OFF" = "#A7A9AC")) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_y_continuous(labels = abs) +
    labs(x = NULL,
         y = "Signaling strength (absolute information flow)",
         fill = "Condition") +
    theme_classic(base_size = 12)

ggsave(file.path(output_dir, "diverging_bar_DMG_tumor_Fig5C.pdf"),
       p_div, width = 4, height = 3.5)

## ══════════════════════════════════════════════════════════════════════════════
## B. LR-pair heatmap – NCAM/NRXN/PTN/CADM (DMG → Neurons) (Fig 5G)
## ══════════════════════════════════════════════════════════════════════════════
# Load the pre-exported net data (exported from CellChat's slot net$df)
# Run once:  write.table(cellchaton@LR$LRsig, "Net_on.txt", sep="\t", quote=FALSE)
#            write.table(cellchatoff@LR$LRsig, "Net_off.txt", sep="\t", quote=FALSE)

df.net.on  <- read.table("data/Net_on.txt",  sep = "\t", header = TRUE,
                          quote = "", stringsAsFactors = FALSE, check.names = FALSE)
df.net.off <- read.table("data/Net_off.txt", sep = "\t", header = TRUE,
                          quote = "", stringsAsFactors = FALSE, check.names = FALSE)

target_celltypes <- c("Neurons", "Excitatory Neurons", "Inhibitory Neurons",
                       "Astrocyte", "Oligodendrocytes")
source_celltype  <- "DMG Tumor"
pathways_hm      <- c("NCAM", "NRXN", "PTN", "CADM")

prepare_hm_matrix <- function(df.net, celltype_targets, src, pws) {
    df.net %>%
        filter(source == src,
               target %in% celltype_targets,
               pathway_name %in% pws) %>%
        group_by(interaction_name_2) %>%
        summarise(prob_sum = sum(prob, na.rm = TRUE), .groups = "drop") %>%
        filter(prob_sum > 0)
}

mat_on  <- prepare_hm_matrix(df.net.on,  target_celltypes, source_celltype, pathways_hm)
mat_off <- prepare_hm_matrix(df.net.off, target_celltypes, source_celltype, pathways_hm)

# Full LR x condition matrix
all_lr <- union(mat_on$interaction_name_2, mat_off$interaction_name_2)
hm_mat <- matrix(0, nrow = length(all_lr), ncol = 2,
                  dimnames = list(all_lr, c("ON", "OFF")))
for (r in mat_on$interaction_name_2)  hm_mat[r, "ON"]  <- mat_on$prob_sum[mat_on$interaction_name_2 == r]
for (r in mat_off$interaction_name_2) hm_mat[r, "OFF"] <- mat_off$prob_sum[mat_off$interaction_name_2 == r]

# Pathway annotation for row coloring
pw_ann_vec <- setNames(
    c(rep("NCAM",4), rep("NRXN",6), rep("PTN",4), rep("CADM",2)),
    all_lr[seq_along(all_lr)]   # simplified; adjust to actual LR names
)

col_fun <- colorRamp2(c(0, max(hm_mat) / 2, max(hm_mat)),
                       c("white", "#FDAE6B", "#D94801"))

ht_lr <- Heatmap(hm_mat,
    name            = "Summed\nprobability",
    col             = col_fun,
    cluster_rows    = FALSE,
    cluster_columns = FALSE,
    row_names_gp    = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    column_title    = "DMG Tumor → Neurons (NCAM/NRXN/PTN/CADM)"
)

pdf(file.path(output_dir, "LR_heatmap_NCAM_NRXN_PTN_CADM_Fig5G.pdf"),
    width = 5, height = 8)
draw(ht_lr)
dev.off()

## ══════════════════════════════════════════════════════════════════════════════
## C. LR-pair dot plot – ON vs OFF (Fig 5H) and OFF-specific (Fig 6F)
## ══════════════════════════════════════════════════════════════════════════════
# Fig 5H: Bubble dot plot comparing LR pair probabilities ON vs OFF
# for NCAM/NRXN/PTN/CADM, DMG sender to Excitatory/Inhibitory Neurons

pathways_fig5H <- c("NCAM", "NRXN", "PTN", "CADM")

p_dot_5H <- netVisual_bubble(
    cellchat,
    sources.use  = "DMG Tumor",
    targets.use  = c("Excitatory Neurons", "Inhibitory Neurons",
                     "Neurons", "Astrocyte", "Oligodendrocytes"),
    signaling    = pathways_fig5H,
    comparison   = c(1, 2),
    angle.x      = 45,
    remove.isolate = TRUE
)
ggsave(file.path(output_dir, "LR_bubble_NCAM_NRXN_PTN_CADM_Fig5H.pdf"),
       p_dot_5H, width = 6, height = 8)

# Fig 6F: OFF-specific pathways bubble plot (Neuron-Neuron)
pathways_fig6F <- c("EPHA", "NT", "OPIOID", "SOMATOSTATIN", "NPY", "EDN")

p_dot_6F <- netVisual_bubble(
    cellchat,
    sources.use  = "Neurons",
    targets.use  = "Neurons",
    signaling    = pathways_fig6F,
    comparison   = c(1, 2),
    angle.x      = 45,
    remove.isolate = TRUE
)
ggsave(file.path(output_dir, "LR_bubble_OFF_pathways_Fig6F.pdf"),
       p_dot_6F, width = 5, height = 7)

cat("\nHeatmap/pathway analysis complete. Outputs saved to:", output_dir, "\n")
