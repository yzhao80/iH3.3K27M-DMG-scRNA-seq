netAnalysis_signalingRole_heatmap2 <- function(
        object,
    signaling = NULL,
    pattern = c("outgoing", "incoming", "all"),
    slot.name = "netP",
    color.use = NULL,
    color.heatmap = "BuGn",
    title = NULL,
    width = 10,
    height = 8,
    font.size = 8,
    font.size.title = 10,
    cluster.rows = FALSE,
    cluster.cols = FALSE,
    top.bar.ylim = NULL,
    right.bar.ylim = NULL
) {
    pattern <- match.arg(pattern)
    
    if (length(slot(object, slot.name)$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
    }
    
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    
    for (i in seq_along(centr)) {
        outgoing[, i] <- centr[[i]]$outdeg
        incoming[, i] <- centr[[i]]$indeg
    }
    
    if (pattern == "outgoing") {
        mat <- t(outgoing)
        legend.name <- "Outgoing"
    } else if (pattern == "incoming") {
        mat <- t(incoming)
        legend.name <- "Incoming"
    } else if (pattern == "all") {
        mat <- t(outgoing + incoming)
        legend.name <- "Overall"
    }
    
    if (is.null(title)) {
        title <- paste0(legend.name, " signaling patterns")
    } else {
        title <- paste0(legend.name, " signaling patterns - ", title)
    }
    
    if (!is.null(signaling)) {
        mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
        mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
        idx <- match(rownames(mat1), signaling)
        mat[idx[!is.na(idx)], ] <- mat1
        dimnames(mat) <- list(signaling, colnames(mat1))
    }
    
    mat.ori <- mat
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    mat[mat == 0] <- NA
    
    ## ---- cell group colors: manual or automatic ----
    group.names <- colnames(mat)
    
    if (is.null(color.use)) {
        color.use <- scPalette(length(group.names))
        names(color.use) <- group.names
    } else {
        # if unnamed vector, assume same order as columns
        if (is.null(names(color.use))) {
            if (length(color.use) != length(group.names)) {
                stop("If `color.use` is unnamed, it must have length equal to the number of cell groups.")
            }
            names(color.use) <- group.names
        } else {
            # if named vector, reorder to match heatmap columns
            if (!all(group.names %in% names(color.use))) {
                missing.groups <- setdiff(group.names, names(color.use))
                stop(
                    "The following cell groups are missing in `color.use`: ",
                    paste(missing.groups, collapse = ", ")
                )
            }
            color.use <- color.use[group.names]
        }
    }
    
    color.heatmap.use <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(n = 9, name = color.heatmap)
    )(100)
    
    df <- data.frame(group = group.names)
    rownames(df) <- group.names
    
    col_annotation <- HeatmapAnnotation(
        df = df,
        col = list(group = color.use),
        which = "column",
        show_legend = FALSE,
        show_annotation_name = FALSE,
        simple_anno_size = grid::unit(0.2, "cm")
    )
    
    ## ---- top barplot ----
    top.bar.values <- colSums(mat.ori)
    if (is.null(top.bar.ylim)) {
        top.bar.ylim <- c(0, max(top.bar.values, na.rm = TRUE))
    }
    
    ha2 <- HeatmapAnnotation(
        Strength = anno_barplot(
            top.bar.values,
            ylim = top.bar.ylim,
            border = FALSE,
            gp = gpar(fill = color.use, col = color.use)
        ),
        show_annotation_name = FALSE
    )
    
    ## ---- right barplot ----
    pSum <- rowSums(mat.ori)
    pSum.original <- pSum
    pSum <- -1 / log(pSum)
    pSum[is.na(pSum)] <- 0
    
    idx1 <- which(is.infinite(pSum) | pSum < 0)
    if (length(idx1) > 0) {
        values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, length.out = length(idx1))
        position <- sort(pSum.original[idx1], index.return = TRUE)$ix
        pSum[idx1] <- values.assign[match(seq_along(idx1), position)]
    }
    
    if (is.null(right.bar.ylim)) {
        right.bar.ylim <- c(0, max(pSum, na.rm = TRUE))
    }
    
    ha1 <- rowAnnotation(
        Strength = anno_barplot(
            pSum,
            ylim = right.bar.ylim,
            border = FALSE
        ),
        show_annotation_name = FALSE
    )
    
    ## ---- legend breaks ----
    if (min(mat, na.rm = TRUE) == max(mat, na.rm = TRUE)) {
        legend.break <- max(mat, na.rm = TRUE)
    } else {
        legend.break <- c(
            round(min(mat, na.rm = TRUE), digits = 1),
            round(max(mat, na.rm = TRUE), digits = 1)
        )
    }
    
    ht1 <- Heatmap(
        mat,
        col = color.heatmap.use,
        na_col = "white",
        name = "Relative strength",
        bottom_annotation = col_annotation,
        top_annotation = ha2,
        right_annotation = ha1,
        cluster_rows = cluster.rows,
        cluster_columns = cluster.cols,
        row_names_side = "left",
        row_names_rot = 0,
        row_names_gp = gpar(fontsize = font.size),
        column_names_gp = gpar(fontsize = font.size),
        width = unit(width, "cm"),
        height = unit(height, "cm"),
        column_title = title,
        column_title_gp = gpar(fontsize = font.size.title),
        column_names_rot = 90,
        heatmap_legend_param = list(
            title_gp = gpar(fontsize = 8, fontface = "plain"),
            title_position = "leftcenter-rot",
            border = NA,
            at = legend.break,
            legend_height = unit(20, "mm"),
            labels_gp = gpar(fontsize = 8),
            grid_width = unit(2, "mm")
        )
    )
    
    return(ht1)
}