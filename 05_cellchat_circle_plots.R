##=============================================================================
## 05_cellchat_circle_plots.R
##
## Project : iH3.3K27M DMG scRNA-seq (Khairkhah et al.)
## Purpose : Generate circle plots showing interaction number and strength
##           for specific signaling pathways (NCAM, NRXN, PTN, CADM) between
##           DMG tumor cells and neuronal populations in ON and OFF conditions.
## Figures : Fig 5D (circle plots per pathway), Sup Fig 7H (differential)
## Input   : cellchaton_analysed.rds, cellchatoff_analysed.rds,
##           cellchat_on_off_merged.rds
##=============================================================================

## ── Libraries ─────────────────────────────────────────────────────────────────
library(CellChat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

## ── Paths ─────────────────────────────────────────────────────────────────────
input_dir  <- "results/04_cellchat/"
output_dir <- "results/05_circle_plots/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ── Load objects ──────────────────────────────────────────────────────────────
cellchaton  <- readRDS(file.path(input_dir, "cellchaton_analysed.rds"))
cellchatoff <- readRDS(file.path(input_dir, "cellchatoff_analysed.rds"))
cellchat    <- readRDS(file.path(input_dir, "cellchat_on_off_merged.rds"))

object.list <- list(ON = cellchaton, OFF = cellchatoff)

## ── Color palette (consistent with UMAP) ─────────────────────────────────────
mycol <- setNames(
    c("#AEC6CF","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
      "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#A6CEE3"),
    c("Astrocyte","DMG Tumor","Endothelial Cells","Ependymal Cells",
      "Epithelial Cells","Excitatory Neurons","Fibroblasts","Immune Cells",
      "Inhibitory Neurons","Microglia","Neurons","Oligodendrocytes",
      "Purkinje Neurons")
)

## ── Helper: circle plot for one pathway, both conditions ─────────────────────
pathway_circle <- function(pathway_name, filename_prefix, width = 12, height = 6) {
    # Compute common edge-weight scale across ON and OFF
    wmax_on  <- max(cellchaton@netP$prob[,,pathway_name],  na.rm = TRUE)
    wmax_off <- max(cellchatoff@netP$prob[,,pathway_name], na.rm = TRUE)
    wmax     <- max(wmax_on, wmax_off, na.rm = TRUE)
    
    png(file.path(output_dir, paste0(filename_prefix, "_circle.png")),
        width = width, height = height, units = "in", res = 600, bg = "white")
    par(mfrow = c(1, 2), xpd = TRUE)
    
    for (nm in names(object.list)) {
        obj <- object.list[[nm]]
        if (pathway_name %in% obj@netP$pathways) {
            netVisual_aggregate(obj,
                signaling      = pathway_name,
                vertex.weight  = as.numeric(table(obj@idents)),
                weight.scale   = TRUE,
                edge.weight.max = wmax,
                color.use      = mycol,
                title.name     = paste(pathway_name, "–", nm))
        } else {
            plot.new()
            title(paste(pathway_name, "–", nm, "(not detected)"))
        }
    }
    dev.off()
    cat("Saved:", file.path(output_dir, paste0(filename_prefix, "_circle.png")), "\n")
}

## ── 1. NCAM pathway (Fig 5D) ─────────────────────────────────────────────────
pathway_circle("NCAM", "NCAM_Fig5D")

## ── 2. NRXN pathway (Fig 5D) ─────────────────────────────────────────────────
pathway_circle("NRXN", "NRXN_Fig5D")

## ── 3. PTN pathway (Fig 5D) ──────────────────────────────────────────────────
pathway_circle("PTN", "PTN_Fig5D")

## ── 4. CADM pathway (Fig 5D) ─────────────────────────────────────────────────
pathway_circle("CADM", "CADM_Fig5D")

## ── 5. OFF-specific pathways (Fig 6): EPHA, NT, OPIOID, SOMATOSTATIN ─────────
off_specific <- c("EPHA", "NT", "OPIOID", "SOMATOSTATIN", "NPY")
for (pw in off_specific) {
    pathway_circle(pw, paste0(pw, "_Fig6"))
}

## ── 6. Aggregate count/strength circles for all cell types (Sup Fig 7H) ──────
# Overall number of interactions
png(file.path(output_dir, "overall_nInteractions_ON_OFF_SupFig7H.png"),
    width = 12, height = 6, units = "in", res = 600, bg = "white")
par(mfrow = c(1, 2), xpd = TRUE)
for (nm in names(object.list)) {
    obj <- object.list[[nm]]
    netVisual_circle(obj@net$count,
        vertex.weight = as.numeric(table(obj@idents)),
        weight.scale  = TRUE,
        label.edge    = FALSE,
        color.use     = mycol,
        title.name    = paste("Number of interactions –", nm))
}
dev.off()

# Overall interaction strength
png(file.path(output_dir, "overall_interactionStrength_ON_OFF_SupFig7H.png"),
    width = 12, height = 6, units = "in", res = 600, bg = "white")
par(mfrow = c(1, 2), xpd = TRUE)
for (nm in names(object.list)) {
    obj <- object.list[[nm]]
    netVisual_circle(obj@net$weight,
        vertex.weight = as.numeric(table(obj@idents)),
        weight.scale  = TRUE,
        label.edge    = FALSE,
        color.use     = mycol,
        title.name    = paste("Interaction strength –", nm))
}
dev.off()

## ── 7. Bar chart: total interactions & strength (Fig 6 inset) ────────────────
# Bar chart showing total number of inferred interactions and total interaction
# strength for ON vs OFF (Fig 6 / Sup Fig 7I)
n_on    <- sum(cellchaton@net$count,  na.rm = TRUE)
n_off   <- sum(cellchatoff@net$count, na.rm = TRUE)
str_on  <- sum(cellchaton@net$weight,  na.rm = TRUE)
str_off <- sum(cellchatoff@net$weight, na.rm = TRUE)

bar_df <- data.frame(
    Condition = rep(c("ON", "OFF"), 2),
    Metric    = rep(c("N interactions", "Interaction strength"), each = 2),
    Value     = c(n_on, n_off, str_on, str_off)
)

p_bar <- ggplot(bar_df, aes(x = Condition, y = Value, fill = Condition)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = round(Value, 3)), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("ON" = "#F1C40F", "OFF" = "#A7A9AC")) +
    facet_wrap(~Metric, scales = "free_y") +
    theme_classic(base_size = 12) +
    labs(x = NULL, y = NULL) +
    NoLegend()

ggsave(file.path(output_dir, "total_interactions_bar_SupFig7I.pdf"),
       p_bar, width = 5, height = 4)

cat("\nCircle plot analysis complete. Outputs saved to:", output_dir, "\n")
