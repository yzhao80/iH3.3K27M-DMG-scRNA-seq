library(Seurat)
library(ggplot2)
library(patchwork)

#----------------------------------------
# 1. format adjusted p-value
#----------------------------------------
fmt_adj_p <- function(p, digits = 2) {
    if (is.na(p)) return("adj.P = NA")
    if (p < 2.2e-16) return("adj.P < 2.2e-16")
    paste0("adj.P = ", formatC(p, format = "e", digits = digits))
}

#----------------------------------------
# 2. make label text for one gene
#----------------------------------------
make_de_label <- function(
        gene,
    de_table,
    p_col = "p_val_adj",
    logfc_col = "avg_log2FC",
    p_digits = 2,
    logfc_digits = 2,
    include_logfc = TRUE
) {
    if (!(gene %in% rownames(de_table))) {
        return("gene not found")
    }
    
    p_adj <- de_table[gene, p_col]
    
    if (include_logfc) {
        logfc <- de_table[gene, logfc_col]
        lab <- paste0(
            fmt_adj_p(p_adj, digits = p_digits),
            "\nlog2FC = ", round(logfc, logfc_digits)
        )
    } else {
        lab <- fmt_adj_p(p_adj, digits = p_digits)
    }
    
    lab
}

#----------------------------------------
# 3. build a list of violin plots
#----------------------------------------
make_violin_plot_list <- function(
        seu,
    genes,
    de_table,
    group.by = "treatment",
    cols = c("ON" = "#F1C40F", "OFF" = "#A7A9AC"),
    pt.size = 0,
    y_lim = c(0, 7),
    label_size = 3.6,
    label_vjust = 1.5,
    p_col = "p_val_adj",
    logfc_col = "avg_log2FC",
    include_logfc = TRUE,
    rotate_x = TRUE
) {
    # keep only genes present in the Seurat object
    genes_present <- genes[genes %in% rownames(seu)]
    
    if (length(genes_present) == 0) {
        stop("None of the requested genes are present in the Seurat object.")
    }
    
    # generate base violin plots
    plist <- VlnPlot(
        seu,
        features = genes_present,
        group.by = group.by,
        pt.size = pt.size,
        combine = FALSE,
        cols = cols
    )
    
    # add labels and formatting
    plist2 <- lapply(seq_along(genes_present), function(i) {
        g <- genes_present[i]
        
        lab <- make_de_label(
            gene = g,
            de_table = de_table,
            p_col = p_col,
            logfc_col = logfc_col,
            include_logfc = include_logfc
        )
        
        p <- plist[[i]] +
            annotate(
                "text",
                x = 1.5,
                y = Inf,
                label = lab,
                vjust = label_vjust,
                size = label_size
            ) +
            coord_cartesian(ylim = y_lim, clip = "off") +
            labs(x = NULL) +
            theme_classic() +
            theme(
                plot.title = element_text(face = "bold", hjust = 0.5),
                plot.margin = margin(t = 22, r = 5, b = 5, l = 5)
            )
        
        if (rotate_x) {
            p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
        
        p
    })
    
    return(plist2)
}

#----------------------------------------
# 4. combine violin plots into one figure
#----------------------------------------
plot_violin_panels <- function(
        seu,
    genes,
    de_table,
    group.by = "treatment",
    cols = c("ON" = "#F1C40F", "OFF" = "#A7A9AC"),
    pt.size = 0,
    y_lim = c(0, 7),
    label_size = 3.6,
    label_vjust = 1.5,
    p_col = "p_val_adj",
    logfc_col = "avg_log2FC",
    include_logfc = TRUE,
    rotate_x = TRUE,
    ncol = 4
) {
    plist2 <- make_violin_plot_list(
        seu = seu,
        genes = genes,
        de_table = de_table,
        group.by = group.by,
        cols = cols,
        pt.size = pt.size,
        y_lim = y_lim,
        label_size = label_size,
        label_vjust = label_vjust,
        p_col = p_col,
        logfc_col = logfc_col,
        include_logfc = include_logfc,
        rotate_x = rotate_x
    )
    
    wrap_plots(plist2, ncol = ncol)
}
