rankNet_patched <- function (object, slot.name = "netP", measure = c("weight", "count"), 
    mode = c("comparison", "single"), comparison = c(1, 2), color.use = NULL, 
    stacked = FALSE, sources.use = NULL, targets.use = NULL, 
    signaling = NULL, pairLR = NULL, signaling.type = NULL, do.stat = FALSE, 
    cutoff.pvalue = 0.05, tol = 0.05, thresh = 0.05, show.raw = FALSE, 
    return.data = FALSE, x.rotation = 90, title = NULL, bar.w = 0.75, 
    font.size = 8, do.flip = TRUE, x.angle = NULL, y.angle = 0, 
    x.hjust = 1, y.hjust = 1, axis.gap = FALSE, ylim = NULL, 
    segments = NULL, tick_width = NULL, rel_heights = c(0.9, 
        0, 0.1)) 
{
    measure <- match.arg(measure)
    mode <- match.arg(mode)
    options(warn = -1)
    object.names <- names(methods::slot(object, slot.name))
    if (measure == "weight") {
        ylabel = "Information flow"
    }
    else if (measure == "count") {
        ylabel = "Number of interactions"
    }
    
    if (mode == "single") {
        stop("This patched version only targets mode = 'comparison'. Use original rankNet() for mode = 'single'.")
    }
    
    ## ------------------------------
    ## mode == "comparison"
    ## ------------------------------
    prob.list <- list()
    pSum <- list()
    pSum.original <- list()
    pair.name <- list()
    idx <- list()
    pSum.original.all <- c()
    object.names.comparison <- c()
    
    for (i in 1:length(comparison)) {
        object.list <- methods::slot(object, slot.name)[[comparison[i]]]
        prob <- object.list$prob
        prob[object.list$pval > thresh] <- 0
        if (measure == "count") {
            prob <- 1 * (prob > 0)
        }
        prob.list[[i]] <- prob
        
        if (!is.null(sources.use)) {
            if (is.character(sources.use)) {
                if (all(sources.use %in% dimnames(prob)[[1]])) {
                    sources.use <- match(sources.use, dimnames(prob)[[1]])
                } else {
                    stop("The input `sources.use` should be cell group names or a numerical vector!")
                }
            }
            idx.t <- setdiff(1:nrow(prob), sources.use)
            prob[idx.t, , ] <- 0
        }
        
        if (!is.null(targets.use)) {
            if (is.character(targets.use)) {
                if (all(targets.use %in% dimnames(prob)[[1]])) {
                    targets.use <- match(targets.use, dimnames(prob)[[2]])
                } else {
                    stop("The input `targets.use` should be cell group names or a numerical vector!")
                }
            }
            idx.t <- setdiff(1:nrow(prob), targets.use)
            prob[, idx.t, ] <- 0
        }
        
        if (sum(prob) == 0) stop("No inferred communications for the input!")
        
        pSum.original[[i]] <- apply(prob, 3, sum)
        if (measure == "weight") {
            pSum[[i]] <- -1/log(pSum.original[[i]])
            pSum[[i]][is.na(pSum[[i]])] <- 0
            idx[[i]] <- which(is.infinite(pSum[[i]]) | pSum[[i]] < 0)
            pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx[[i]]])
        } else {
            pSum[[i]] <- pSum.original[[i]]
        }
        
        pair.name[[i]] <- names(pSum.original[[i]])
        object.names.comparison <- c(object.names.comparison, object.names[comparison[i]])
    }
    
    if (measure == "weight") {
        values.assign <- seq(max(unlist(pSum)) * 1.1, max(unlist(pSum)) * 1.5, length.out = length(unlist(idx)))
        position <- sort(pSum.original.all, index.return = TRUE)$ix
        for (i in 1:length(comparison)) {
            if (i == 1) {
                pSum[[i]][idx[[i]]] <- values.assign[match(1:length(idx[[i]]), position)]
            } else {
                pSum[[i]][idx[[i]]] <- values.assign[match(length(unlist(idx[1:i - 1])) + 1:length(unlist(idx[1:i])), position)]
            }
        }
    }
    
    pair.name.all <- as.character(unique(unlist(pair.name)))
    df <- list()
    for (i in 1:length(comparison)) {
        df[[i]] <- data.frame(name = pair.name.all, contribution = 0, contribution.scaled = 0,
            group = object.names[comparison[i]], row.names = pair.name.all)
        df[[i]][pair.name[[i]], 3] <- pSum[[i]]
        df[[i]][pair.name[[i]], 2] <- pSum.original[[i]]
    }
    
    contribution.relative <- list()
    for (i in 1:(length(comparison) - 1)) {
        contribution.relative[[i]] <- as.numeric(format(df[[length(comparison) - i + 1]]$contribution/df[[1]]$contribution, digits = 1))
        contribution.relative[[i]][is.na(contribution.relative[[i]])] <- 0
    }
    names(contribution.relative) <- paste0("contribution.relative.", 1:length(contribution.relative))
    
    for (i in 1:length(comparison)) {
        for (j in 1:length(contribution.relative)) {
            df[[i]][[names(contribution.relative)[j]]] <- contribution.relative[[j]]
        }
    }
    
    df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution
    if (length(comparison) == 2) {
        idx <- with(df[[1]], order(-contribution.relative.1, contribution, -contribution.data2))
    } else if (length(comparison) == 3) {
        idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2, contribution, -contribution.data2))
    } else if (length(comparison) == 4) {
        idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2, -contribution.relative.3, contribution, -contribution.data2))
    } else {
        idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2, -contribution.relative.3, -contribution.relative.4, contribution, -contribution.data2))
    }
    
    for (i in 1:length(comparison)) {
        df[[i]] <- df[[i]][idx, ]
        df[[i]]$name <- factor(df[[i]]$name, levels = as.character(df[[i]]$name))
    }
    df[[1]]$contribution.data2 <- NULL
    df <- do.call(rbind, df)
    df$group <- factor(df$group, levels = rev(object.names.comparison))
    
    ## --- CUSTOM BAR COLORS (your mapping) ---
    ## We lock in your colors for 2-group comparison unless user explicitly passes color.use
    if (is.null(color.use)) {
        if (length(comparison) == 2) {
            color.use <- c("#A7A9AC", "#F1C40F")  # "red" -> #F1C40F, "green" -> #A7A9AC
        } else {
            color.use <- ggPalette(length(comparison))
        }
    }
    
    ## --- STAT TEST (same as original; fixes a small bug: j -> i not needed) ---
    if (do.stat && length(comparison) == 2) {
        df$pvalues <- 0
        for (k in 1:length(pair.name.all)) {
            if (nrow(prob.list[[1]]) != nrow(prob.list[[2]])) {
                stop("Statistical test is not applicable to datasets with different cellular compositions! Please set `do.stat = FALSE`")
            }
            prob.values <- matrix(0, nrow = nrow(prob.list[[1]]) * nrow(prob.list[[1]]), ncol = length(comparison))
            for (j in 1:length(comparison)) {
                if (pair.name.all[k] %in% pair.name[[j]]) {
                    prob.values[, j] <- as.vector(prob.list[[j]][, , pair.name.all[k]])
                } else {
                    prob.values[, j] <- NA
                }
            }
            prob.values <- prob.values[rowSums(prob.values, na.rm = TRUE) != 0, , drop = FALSE]
            if (nrow(prob.values) > 3 && sum(is.na(prob.values)) == 0) {
                pvalues <- wilcox.test(prob.values[, 1], prob.values[, 2], paired = TRUE)$p.value
            } else {
                pvalues <- 0
            }
            pvalues[is.na(pvalues)] <- 0
            df$pvalues[df$name == pair.name.all[k]] <- pvalues
        }
    }
    
    ## --- REMOVE "NON-SIGNIFICANT" PATHWAYS (the ones that would have been black) ---
    if (length(comparison) == 2) {
        df_one <- df[!duplicated(df$name), ]  # one row per pathway (same relative/pvalue repeated)
        if (do.stat) {
            keep <- ((df_one$contribution.relative.1 < (1 - tol)) | (df_one$contribution.relative.1 > (1 + tol))) &
                (df_one$pvalues < cutoff.pvalue)
        } else {
            keep <- ((df_one$contribution.relative.1 < (1 - tol)) | (df_one$contribution.relative.1 > (1 + tol)))
        }
        keep_names <- as.character(df_one$name[keep])
        df <- df[df$name %in% keep_names, , drop = FALSE]
        df$name <- factor(df$name, levels = levels(df_one$name)[levels(df_one$name) %in% keep_names])
    }
    
    ## Remove zero-sum pathways (original behavior)
    for (i in 1:length(pair.name.all)) {
        df.t <- df[df$name == pair.name.all[i], "contribution"]
        if (length(df.t) > 0 && sum(df.t) == 0) {
            df <- df[-which(df$name == pair.name.all[i]), ]
        }
    }
    
    ## Optional filtering by signaling / pairLR (original behavior)
    if ((slot.name == "netP") && (!is.null(signaling))) {
        df <- subset(df, name %in% signaling)
    } else if ((slot.name == "netP") && (!is.null(pairLR))) {
        stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
        df <- subset(df, name %in% pairLR)
    } else if ((slot.name == "net") && (!is.null(signaling))) {
        stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }
    
    ## --- ALL PATHWAY NAMES BLACK ---
    axis_text_colors <- rep("black", length(levels(df$name)))
    
    ## Plot
    if (stacked) {
        gg <- ggplot(df, aes(x = name, y = contribution, fill = group)) +
            geom_bar(stat = "identity", width = bar.w, position = "fill")
        if (measure == "weight") {
            gg <- gg + xlab("") + ylab("Relative information flow")
        } else {
            gg <- gg + xlab("") + ylab("Relative number of interactions")
        }
        gg <- gg + geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", size = 0.5)
    } else {
        if (show.raw) {
            gg <- ggplot(df, aes(x = name, y = contribution, fill = group)) +
                geom_bar(stat = "identity", width = bar.w, position = position_dodge(0.8)) +
                xlab("") + ylab(ylabel)
        } else {
            gg <- ggplot(df, aes(x = name, y = contribution.scaled, fill = group)) +
                geom_bar(stat = "identity", width = bar.w, position = position_dodge(0.8)) +
                xlab("") + ylab(ylabel)
        }
        if (axis.gap) {
            gg <- gg + theme_bw() + theme(panel.grid = element_blank())
            gg.gap::gg.gap(gg, ylim = ylim, segments = segments, tick_width = tick_width, rel_heights = rel_heights)
        }
    }
    
    gg <- gg + CellChat_theme_opts() + theme_classic()
    
    if (do.flip) {
        gg <- gg + coord_flip() +
            theme(axis.text.y = element_text(colour = axis_text_colors))
        if (is.null(x.angle)) x.angle = 0
    } else {
        if (is.null(x.angle)) x.angle = 45
        gg <- gg + scale_x_discrete(limits = rev) +
            theme(axis.text.x = element_text(colour = rev(axis_text_colors)))
    }
    
    gg <- gg + theme(axis.text = element_text(size = font.size),
        axis.title.y = element_text(size = font.size))
    
    gg <- gg + scale_fill_manual(name = "", values = color.use) +
        guides(fill = guide_legend(reverse = T)) +
        theme(axis.text.x = element_text(angle = x.angle, hjust = x.hjust),
            axis.text.y = element_text(angle = y.angle, hjust = y.hjust))
    
    if (!is.null(title)) {
        gg <- gg + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    }
    
    if (return.data) {
        df$contribution <- abs(df$contribution)
        df$contribution.scaled <- abs(df$contribution.scaled)
        return(list(signaling.contribution = df, gg.obj = gg))
    } else {
        return(gg)
    }
}
