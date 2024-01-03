
library(ggrepel)

cluster1.markers <- markers %>% 
  dplyr::filter(.,cluster == 'H') %>% 
  mutate(Difference = pct.1 - pct.2)

log2FC = 1
padj = 0.05 

cluster1.markers$threshold="ns";
cluster1.markers[which(cluster1.markers$avg_log2FC  > log2FC & cluster1.markers$p_val_adj <padj),]$threshold="up";
cluster1.markers[which(cluster1.markers$avg_log2FC  < (-log2FC) & cluster1.markers$p_val_adj < padj),]$threshold="down";
cluster1.markers$threshold=factor(cluster1.markers$threshold, levels=c('down','ns','up'))

ggplot(cluster1.markers, aes(x=Difference, y= -log10(p_val_adj) , color = threshold)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=c( "blue","grey","red") ) + 
  geom_label_repel(data=subset(cluster1.markers, avg_log2FC >= 1 & Difference >= 0.2 & p_val_adj <= 0.05), 
                   aes(label=gene),  #添加label
                   color="black", #设置label中标签的颜色
                   segment.colour = "black",#设置label框的颜色
                   label.padding = 0.1, 
                   #max.overlaps = 200,
                   segment.size = 0.3,  #框的大小
                   size=4)+
  geom_label_repel(data=subset(cluster1.markers, avg_log2FC <= -1 & Difference <= -0.2 & p_val_adj <= 0.05), 
                   aes(label=gene), label.padding = 0.1, 
                   color="black",
                   segment.colour = "black",
                   segment.size = 0.3, size=4)+
  geom_vline(xintercept = 0.0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  theme_classic()

VolcanoPlot <- function(srt, group_by = NULL, test.use = "wilcox", DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
                        x_metric = "diff_pct", palette = "RdBu", palcolor = NULL, pt.size = 1, pt.alpha = 1,
                        cols.highlight = "black", sizes.highlight = 1, alpha.highlight = 1, stroke.highlight = 0.5,
                        nlabel = 5, features_label = NULL, label.fg = "black", label.bg = "white", label.bg.r = 0.1, label.size = 4,
                        aspect.ratio = NULL, xlab = x_metric, ylab = "-log10(p-adjust)",
                        theme_use = "theme_scp", theme_args = list(),
                        combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE) {
  if (is.null(group_by)) {
    group_by <- "custom"
  }
  slot <- paste0("DEtest_", group_by)
  if (!slot %in% names(srt@tools) || length(grep(pattern = "AllMarkers", names(srt@tools[[slot]]))) == 0) {
    stop("Cannot find the DEtest result for the group '", group_by, "'. You may perform RunDEtest first.")
  }
  index <- grep(pattern = paste0("AllMarkers_", test.use), names(srt@tools[[slot]]))[1]
  if (is.na(index)) {
    stop("Cannot find the 'AllMarkers_", test.use, "' in the DEtest result.")
  }
  de <- names(srt@tools[[slot]])[index]
  de_df <- srt@tools[[slot]][[de]]
  de_df[, "diff_pct"] <- de_df[, "pct.1"] - de_df[, "pct.2"]
  de_df[, "-log10padj"] <- -log10(de_df[, "p_val_adj"])
  de_df[, "DE"] <- FALSE
  de_df[with(de_df, eval(rlang::parse_expr(DE_threshold))), "DE"] <- TRUE
  
  x_upper <- quantile(de_df[["avg_log2FC"]][is.finite(de_df[["avg_log2FC"]])], c(0.99, 1))
  x_lower <- quantile(de_df[["avg_log2FC"]][is.finite(de_df[["avg_log2FC"]])], c(0.01, 0))
  x_upper <- ifelse(x_upper[1] > 0, x_upper[1], x_upper[2])
  x_lower <- ifelse(x_lower[1] < 0, x_lower[1], x_lower[2])
  if (x_upper > 0 & x_lower < 0) {
    value_range <- min(abs(c(x_upper, x_lower)), na.rm = TRUE)
    x_upper <- value_range
    x_lower <- -value_range
  }
  
  de_df[, "border"] <- FALSE
  de_df[de_df[["avg_log2FC"]] > x_upper, "border"] <- TRUE
  de_df[de_df[["avg_log2FC"]] > x_upper, "avg_log2FC"] <- x_upper
  de_df[de_df[["avg_log2FC"]] < x_lower, "border"] <- TRUE
  de_df[de_df[["avg_log2FC"]] < x_lower, "avg_log2FC"] <- x_lower
  
  de_df[, "y"] <- -log10(de_df[, "p_val_adj"])
  if (x_metric == "diff_pct") {
    de_df[, "x"] <- de_df[, "diff_pct"]
    de_df[de_df[, "avg_log2FC"] < 0, "y"] <- -de_df[de_df[, "avg_log2FC"] < 0, "y"]
    de_df <- de_df[order(abs(de_df[, "avg_log2FC"]), decreasing = FALSE, na.last = FALSE), , drop = FALSE]
  } else if (x_metric == "avg_log2FC") {
    de_df[, "x"] <- de_df[, "avg_log2FC"]
    de_df[de_df[, "diff_pct"] < 0, "y"] <- -de_df[de_df[, "diff_pct"] < 0, "y"]
    de_df <- de_df[order(abs(de_df[, "diff_pct"]), decreasing = FALSE, na.last = FALSE), , drop = FALSE]
  }
  de_df[, "distance"] <- de_df[, "x"]^2 + de_df[, "y"]^2
  
  plist <- list()
  for (group in levels(de_df[["group1"]])) {
    df <- de_df[de_df[["group1"]] == group, , drop = FALSE]
    if (nrow(df) == 0) {
      next
    }
    x_nudge <- diff(range(df$x)) * 0.05
    df[, "label"] <- FALSE
    if (is.null(features_label)) {
      df[df[["y"]] >= 0, ][head(order(df[df[["y"]] >= 0, "distance"], decreasing = TRUE), nlabel), "label"] <- TRUE
      df[df[["y"]] < 0, ][head(order(df[df[["y"]] < 0, "distance"], decreasing = TRUE), nlabel), "label"] <- TRUE
    } else {
      df[df[["gene"]] %in% features_label, "label"] <- TRUE
    }
    jitter <- position_jitter(width = 0.2, height = 0.2, seed = 11)
    color_by <- ifelse(x_metric == "diff_pct", "avg_log2FC", "diff_pct")
    p <- ggplot() +
      geom_point(data = df[!df[["DE"]] & !df[["border"]], , drop = FALSE], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha) +
      geom_point(data = df[!df[["DE"]] & df[["border"]], , drop = FALSE], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha, position = jitter) +
      geom_point(data = df[df[["DE"]] & !df[["border"]], , drop = FALSE], aes(x = x, y = y), color = cols.highlight, size = sizes.highlight + stroke.highlight, alpha = alpha.highlight) +
      geom_point(data = df[df[["DE"]] & df[["border"]], , drop = FALSE], aes(x = x, y = y), color = cols.highlight, size = sizes.highlight + stroke.highlight, alpha = alpha.highlight, position = jitter) +
      geom_point(data = df[df[["DE"]] & !df[["border"]], , drop = FALSE], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha) +
      geom_point(data = df[df[["DE"]] & df[["border"]], , drop = FALSE], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha, position = jitter) +
      geom_hline(yintercept = 0, color = "black", linetype = 1) +
      geom_vline(xintercept = 0, color = "grey", linetype = 2) +
      geom_text_repel(
        data = df[df[["label"]], , drop = FALSE], aes(x = x, y = y, label = gene),
        min.segment.length = 0, max.overlaps = 100, segment.colour = "grey40",
        color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, force = 20,
        nudge_x = ifelse(df[df[["label"]], "y"] >= 0, -x_nudge, x_nudge)
      ) +
      labs(x = xlab, y = ylab) +
      scale_color_gradientn(
        name = ifelse(x_metric == "diff_pct", "log2FC", "diff_pct"), colors = palette_scp(palette = palette, palcolor = palcolor),
        values = rescale(unique(c(min(c(df[, color_by], 0), na.rm = TRUE), 0, max(df[, color_by], na.rm = TRUE)))),
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)
      ) +
      scale_y_continuous(labels = abs) +
      facet_wrap(~group1) +
      do.call(theme_use, theme_args) +
      theme(aspect.ratio = aspect.ratio)
    plist[[group]] <- p
  }
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}
