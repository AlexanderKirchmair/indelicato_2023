
## R functions for plotting and visualization

ggvolcano <- function(data, x = NULL, y = NULL, color = NULL, label = NULL, shape = NULL,
                      nlabels = NULL, lab_size = 12, lab_ix = NULL, repel = 1.5, attract = NULL, box.padding = 0.5, max_overlaps = Inf, seed = 123,
                      ptres = 0.05, clip = FALSE, symlim = TRUE, expand = c(0,0), nbreaks_x = 7, nbreaks_y = 7,
                      xlim = NULL, ylim = NULL,
                      color_up = "#eb9d0e", color_down = "#146bc7", color_nonsig = "#4d4d4d",
                      title = NULL, title_size = NULL, point_size = 2, scale_size = FALSE, axis_size = NULL, leg_size = NULL,
                      lwd = 0.8, at_zero = FALSE, ...){
  
  ### Function to plot volcano plots using ggplot.
  
  data <- as.data.frame(data)
  
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  
  if (rlang::quo_is_null(x)){ x <- rlang::sym(grep("lfc|log2FoldChange|logFC|log2FC|nes", names(data), value = TRUE, ignore.case = TRUE)[1]) }
  if (rlang::quo_is_null(y)){ y <- rlang::sym(grep("padj|fdr", names(data), value = TRUE, ignore.case = TRUE)[1]) }
  
  data$x <- data[[rlang::as_name(x)]]
  data$y <- -log10(data[[rlang::as_name(y)]])
  data <- data[!is.na(data$x) & !is.na(data$y),]
  
  data$xtmp <- data$x
  data$xtmp[is.infinite(data$xtmp)] <- max(abs(data$xtmp[!is.infinite(data$xtmp)])) * sign(data$xtmp[is.infinite(data$xtmp)])
  data$ytmp <- data$y
  data$ytmp[is.infinite(data$ytmp)] <- max(abs(data$ytmp[!is.infinite(data$ytmp)])) * sign(data$ytmp[is.infinite(data$ytmp)])
  data$score <- abs(as.numeric(scale(data$xtmp, center = FALSE))) + abs(as.numeric(scale(data$ytmp, center = FALSE)))
  data$score[is.na(data$score)] <- 0
  
  data$class <- "not signif."
  data$class[data[[rlang::as_name(y)]] <= ptres & data$x > 0] <- "up"
  data$class[data[[rlang::as_name(y)]] <= ptres & data$x < 0] <- "down"
  
  data$score[data$class == "not signif."] <- data$score[data$class == "not signif."] * 0.001
  
  if (is.null(title_size)) title_size <- lab_size
  if (is.null(axis_size)) axis_size <- lab_size
  if (is.null(leg_size)) leg_size <- lab_size
  
  
  shape <- rlang::enquo(shape)
  
  
  # LABELS
  
  label <- rlang::enquo(label)
  
  if (rlang::quo_is_null(label)){
    data[["label"]] <- rownames(data)
  } else {
    data[["label"]] <- data[[rlang::as_name(label)]]
  }
  
  if (is.null(lab_ix)){
    
    data <- data[order(data$score, decreasing = TRUE),]
    
    if (is.null(nlabels)){ nlabels <- min(20, ceiling(nrow(data)/10)) }
    if (is.infinite(nlabels)){ nlabels <- nrow(data) }
    data$do_label <- FALSE
    nlabels_left <- nlabels_right <- 0
    if (nrow(subset(data, x < 0)) > 0) nlabels_left <- ceiling(nlabels/2 * max(subset(data, x < 0)$score, na.rm = TRUE) / max(data$score, na.rm = TRUE))
    if (nrow(subset(data, x > 0)) > 0) nlabels_right <- ceiling(nlabels/2 * max(subset(data, x > 0)$score, na.rm = TRUE) / max(data$score, na.rm = TRUE))
    if (is.na(nlabels_left)) nlabels_left <- 0
    if (is.na(nlabels_right)) nlabels_right <- 0
    
    data$do_label[data$x < 0][1:nlabels_left] <- TRUE
    data$do_label[data$x > 0][1:nlabels_right] <- TRUE
    data$do_label[is.na(data$do_label)] <- FALSE
    
    if (sum(data$do_label) < nlabels){ data$do_label[!data$do_label][1:(nlabels-sum(data$do_label))] <- TRUE }
    data$label[!data$do_label] <- ""
    data$do_label[data$label == ""] <- FALSE
    
    data$do_label[data$class == "not signif."] <- FALSE
    
  } else {
    
    data$do_label <- FALSE
    data[lab_ix,]$do_label <- TRUE
    
  }
  
  
  
  # COLORS
  
  color <- rlang::enquo(color)
  
  if (rlang::quo_is_null(color)){
    color <- rlang::sym("class")
    colorvals <-  c("up" = color_up, "down" = color_down, "not signif." = color_nonsig)
  } else {
    col_levels <- unique(data[[rlang::as_name(color)]])
    colorvals <-  setNames(scales::muted(rainbow(length(col_levels))), col_levels)
  }
  
  
  # LIMITS
  
  sigdata <- subset(data, class != "not signif.")
  xylimits <- list(xlim = getLimits(sigdata$xtmp, clip = clip, expand = expand[1]), ylim = getLimits(sigdata$ytmp, clip = clip, expand = expand[2], negative = FALSE))
  if (symlim == TRUE){ xylimits$xlim <- c("min" = -max(abs(xylimits$xlim)), "max" = max(abs(xylimits$xlim))) }
  data$xorg <- data$x
  data$yorg <- data$y
  
  if (!is.null(xlim)) xylimits$xlim <- setNames(xlim, c("min", "max"))
  if (!is.null(ylim)) xylimits$ylim <- setNames(ylim, c("min", "max"))
  
  # X breaks
  xclip_min <- any(naf(data$xorg < xylimits$xlim["min"]))
  xclip_max <- any(naf(data$xorg > xylimits$xlim["max"]))
  
  xbreaks <- scales::pretty_breaks(n = nbreaks_x)(xylimits$xlim, n = nbreaks_x)
  if (clip & xclip_min){
    xylimits$xlim["min"] <- min(xbreaks)
    xclip_min <- any(data$xorg < xylimits$xlim["min"])
  }
  if (clip & xclip_max){
    xylimits$xlim["max"] <- max(xbreaks)
    xclip_max <- any(data$xorg > xylimits$xlim["max"])
  }
  
  
  xylimits$xlim <- xylimits$xlim + c(-diff(xylimits$xlim), diff(xylimits$xlim)) * c(!xclip_min, !xclip_max)*0.02
  
  names(xbreaks) <- as.character(xbreaks)
  if (xclip_min){ names(xbreaks)[1] <- paste0("<", xbreaks[1]) }
  if (xclip_max){ names(xbreaks)[length(xbreaks)] <- paste0(">", xbreaks[length(xbreaks)]) }
  
  
  # Y breaks
  yclip_min <- any(naf(data$yorg < xylimits$ylim["min"]))
  yclip_max <- any(naf(data$yorg > xylimits$ylim["max"]))
  
  ybreaks <- scales::pretty_breaks(n = nbreaks_y)(xylimits$ylim, n = nbreaks_y)
  
  if (clip & yclip_min){
    xylimits$ylim["min"] <- min(ybreaks)
    yclip_min <- any(data$xorg < xylimits$ylim["min"])
  }
  if (clip & yclip_max){
    xylimits$ylim["max"] <- max(ybreaks)
    yclip_max <- any(data$xorg > xylimits$ylim["max"])
  }
  
  
  xylimits$ylim <- xylimits$ylim + c(-diff(xylimits$ylim), diff(xylimits$ylim)) * c(!yclip_min & !at_zero, !yclip_max)*c(0.01, 0.05)
  
  names(ybreaks) <- as.character(ybreaks)
  if (yclip_min){ names(ybreaks)[1] <- paste0("<", ybreaks[1]) }
  if (yclip_max){ names(ybreaks)[length(ybreaks)] <- paste0(">", ybreaks[length(ybreaks)]) }
  
  
  data$x[ data$x < xylimits$xlim["min"] ] <- xylimits$xlim["min"]
  data$x[ data$x > xylimits$xlim["max"] ] <- xylimits$xlim["max"]
  data$y[ data$y < xylimits$ylim["min"] ] <- xylimits$ylim["min"]
  data$y[ data$y > xylimits$ylim["max"] ] <- xylimits$ylim["max"]
  
  
  
  ### GGPLOT ###
  
  data <- data[order(data$score, decreasing = FALSE),]
  
  gg <- data %>% ggplot2::ggplot(ggplot2::aes(x = x, y = y, label = label, color = !!color, shape = !!shape, ...))
  
  gg %<>% + ggplot2::theme_bw(base_size = 20)
  gg %<>% + ggplot2::theme(text = ggplot2::element_text(color = "black", size = lab_size),
                           rect = ggplot2::element_rect(color = "black", size = lwd),
                           line = ggplot2::element_line(size = lwd),
                           legend.text = ggplot2::element_text(color = "black", size = leg_size),
                           legend.title = ggplot2::element_text(color = "black", size = leg_size),
                           panel.grid.minor = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_line(size = lwd, color = rgb(0.9,0.9,0.9)),
                           panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = lwd),
                           strip.background = ggplot2::element_blank(),
                           strip.text = ggplot2::element_text(color = "black", size = title_size),
                           axis.ticks = ggplot2::element_line(color = "black", size = lwd),
                           axis.line = ggplot2::element_blank(),
                           plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
                           plot.title = ggplot2::element_text(size = title_size, hjust = 0.5, lineheight = 1.5),
                           axis.title = ggplot2::element_text(size = axis_size, face = "bold"),
                           axis.text = ggplot2::element_text(size = axis_size, color = "black"))
  
  if (!is.null(ptres)){ gg %<>% + ggplot2::geom_hline(yintercept = -log10(ptres), linetype = "dashed", color = rgb(0.3,0.3,0.3)) }
  
  # points
  if (scale_size == FALSE){
    gg %<>% + ggplot2::geom_point(size = point_size, alpha = 0.8)
  } else {
    gg %<>% + ggplot2::geom_point(aes(size = score), alpha = 0.8)
    gg %<>% + ggplot2::scale_size_continuous(range = c(point_size/5, point_size*2), guide = "none")
  }
  
  
  
  gg %<>% + ggplot2::scale_colour_manual(values = colorvals)
  gg %<>% + ggplot2::labs(title = title, y = paste0("-log10 ", rlang::as_name(y)), x = rlang::as_name(x), size = "none")
  
  
  gg %<>% + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0,0)),
                                        limits = xylimits$xlim,
                                        breaks = xbreaks,
                                        labels = names(xbreaks))
  
  gg %<>% + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0)),
                                        limits = xylimits$ylim,
                                        breaks = ybreaks,
                                        labels = names(ybreaks))
  
  # point labels
  if (is.null(attract)) attract <- sqrt(repel)
  gg %<>% + ggrepel::geom_text_repel(data = subset(data, do_label == TRUE),
                                     size = lab_size/ggplot2:::.pt,
                                     seed = seed,
                                     xlim = xylimits$xlim - c(-diff(xylimits$xlim), diff(xylimits$xlim))*0.18,
                                     ylim = xylimits$ylim - c(-diff(xylimits$ylim)*0.3, diff(xylimits$ylim)*0.02),
                                     force = repel, force_pull = attract,  max.overlaps = max_overlaps,
                                     point.padding = 0.35, box.padding = box.padding, max.time = 30, max.iter = 10^6, min.segment.length = 0, vjust = 0, color = rgb(0.0,0.0,0.0), segment.alpha = 0.6)
  
  gg %<>% + ggplot2::coord_cartesian(clip = "off")
  
  return(gg)
}




getLimits <- function(x, clip = TRUE, expand = 1, negative = TRUE){
  
  x <- x[!is.na(x)]
  x <- x + x*expand
  
  if (clip == TRUE){
    h <- hist(x, plot = FALSE, breaks = 30)
    xd <- h$counts > 3
    xmin <- h$breaks[which(xd)[1]]
    xmax <- rev(h$breaks)[which(rev(xd))[1]]
    
  } else {
    xmin <- NA
    xmax <- NA
  }
  
  if (is.na(xmax)) xmax <- max(x) %>% roundup(., roundup(-log10(abs(.)))) # upper
  if (is.na(xmin)) xmin <- min(x) %>% rounddown(., roundup(-log10(abs(.)))) # lower
  
  if (is.na(xmin)){xmin <- -0.1 * xmax}
  if (is.na(xmax)){xmax <- -0.1 * xmin}
  
  if (is.na(xmax) & is.na(xmin)){
    xmin <- -1
    xmax <- 1
  }
  
  
  res <- c("min" = xmin, "max" = xmax)
  if (negative == FALSE) res[res < 0] <- 0
  
  res
}



ggpca <- function(data, design = NULL, mapping =  ggplot2::aes(), center = TRUE, scale = TRUE, na = NULL, n = NULL, label = FALSE, digits = 2, colors = NULL, size = 3.5, title = NULL, return_data = FALSE, ...){
  
  # annotation data
  if (is.null(design)) design <- data.frame(row.names = rownames(data))
  
  # missing values
  if (any(is.na(data))){
    if (na == "impute"){
      ncp <-  missMDA::estim_ncpPCA(data, ncp.min = 2, ncp.max = max(round(ncol(data)/3), 3))
      data <- missMDA::imputePCA(data, ncp = ncp$ncp)$completeObs
    } else if (na == "omit"){
      data <- na.omit(data)
    } else if (is.numeric(na)){
      data[is.na(data)] <- na
    } else {
      stop("Please select a method for dealing with missing values!")
    }
  }
  
  # select top-variance features
  if (!is.null(n)){
    if (n <= 1) n <- round(nrow(data) * n)
    vars <- apply(data, 1, var)
    data <- data[order(abs(vars), decreasing = TRUE)[1:n],]
  }
  
  # remove zero-variance features
  vars <- apply(data, 1, var, na.rm = TRUE)
  data <- data[abs(vars) > .Machine$double.eps & !is.na(vars),]
    
  # run pca
  pcares <- stats::prcomp(t(data), center = center, scale. = scale, ...)
  if (return_data == TRUE) return(pcares)
  
  # plotting data
  pca_df <- data.frame(pcares$x, design[rownames(pcares$x),,drop = FALSE])
  variance_explained <- summary(pcares)$importance["Proportion of Variance", ]
  pca_df$label <- rownames(pca_df)
  
  # ggplot
  base_aes <- ggplot2::aes(x = PC1, y = PC2, label = label)
  base_aes[names(mapping)] <- mapping
  
  gg <- ggplot2::ggplot(data = pca_df, mapping = base_aes) +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::theme(panel.border =  ggplot2::element_rect(color = "black", fill = NA, size = 1),
                   axis.text =  ggplot2::element_text(color = "black"),
                   axis.ticks =  ggplot2::element_line(color = "black")) +
    ggplot2::geom_point(size = size) +
    ggplot2::xlab(paste0("PC1 (", round(100*variance_explained[1], digits), "%)")) +
    ggplot2::ylab(paste0("PC2 (", round(100*variance_explained[2], digits), "%)")) +
    ggplot2::scale_x_continuous(expand =  ggplot2::expansion(mult = c(0.2,0.2))) +
    ggplot2::scale_y_continuous(expand =  ggplot2::expansion(mult = c(0.2,0.2)))
  
  if (label == TRUE) gg <- gg +  ggplot2::geom_text_repel(size = 5, show.legend = F, min.segment.length = 2)
  if (!is.null(colors)) gg <- gg +  ggplot2::scale_color_manual(values = colors[[rlang::as_name(base_aes[["colour"]])]])
  if (!is.null(title)) gg <- gg + ggplot2::ggtitle(title)
  
  gg
  
}



getColors <- function(dataframe) 
{
  if (class(dataframe) != "data.frame") 
    dataframe <- data.frame(dataframe)
  if (is.null(dataframe)) 
    return(NULL)
  if (ncol(dataframe) == 0) 
    return(NULL)
  colors <- list()
  dataframe.discrete <- dataframe[, sapply(dataframe, is.factor) | 
                                    sapply(dataframe, is.character), drop = FALSE]
  n <- sapply(dataframe.discrete, function(x) length(unique(x)))
  colors <- genPalettes(n = NULL, length_each = n)
  colors <- lapply(seq_along(colors), function(i) {
    tmp <- colors[[i]]
    names(tmp) <- unique(dataframe.discrete[[i]])
    tmp
  })
  names(colors) <- colnames(dataframe.discrete)
  colors
}


genPalettes <- function (n = 1, length_each = 3, dist = 0.75, saturation = c(0.5, 0.7), lightness = c(0.5, 0.7), cvd = "protan", cvd_severity = 0.2, seed = 123, ...){
  
  if (is.null(n)) {
    n <- length(length_each)
  }
  if (length(length_each) == 1 & n != 1) {
    length_each <- rep(length_each, n)
  }
  if (length(length_each) != n) {
    stop("Please provide the number of colors to generate in each palette!")
  }
  rel <- (0:n)/(n)
  drel <- mean(diff(rel))
  sep <- 1 - dist
  init <- scales::rescale(rel, from = c(min(rel) - drel * sep, 
                                        max(rel) + drel * sep), to = c(1, 359))
  d <- mean(diff(init))
  set.seed(seed)
  init <- init[-length(init)] + d * runif(1)
  spaces <- lapply(init, function(tmp) {
    list(h = c(tmp - d * sep, tmp + d * sep), s = saturation, 
         l = lightness)
  })
  pals <- lapply(1:n, function(i) {
    qualpalr::qualpal(n = length_each[i], colorspace = spaces[[i]], 
                      cvd_severity = cvd_severity, cvd = cvd, ...)$hex
  })
  pals <- lapply(pals, function(pal) {
    ix <- order(grDevices::rgb2hsv(grDevices::col2rgb(pal))["v", 
    ], decreasing = TRUE)
    pal[ix]
  })
  names(pals) <- names(length_each)
  pals
}



theme_basic <- function (base_size = 18, base_family = "", base_line_size = base_size/22, 
          base_rect_size = base_size/22, base_color = "black", grid = FALSE, 
          ...){
  
  th1 <- ggplot2::theme_bw(base_size = base_size, base_family = base_family, 
                           base_line_size = base_line_size, base_rect_size = base_rect_size)
  th2 <- ggplot2::theme(line = ggplot2::element_line(colour = base_color, 
                                                     size = base_line_size), rect = ggplot2::element_rect(colour = base_color, 
                                                                                                          fill = NA, size = base_rect_size), text = ggplot2::element_text(colour = base_color, 
                                                                                                                                                                          size = base_size, family = base_family), title = ggplot2::element_text(colour = base_color), 
                        axis.text = ggplot2::element_text(colour = base_color, 
                                                          size = ggplot2::rel(0.75)), axis.ticks = ggplot2::element_line(colour = base_color), 
                        axis.line = ggplot2::element_line(colour = base_color), 
                        legend.text = ggplot2::element_text(colour = base_color), 
                        legend.title = ggplot2::element_text(colour = base_color), 
                        legend.title.align = 0, panel.border = ggplot2::element_blank(), 
                        panel.grid = ggplot2::element_blank(), plot.title = ggplot2::element_text(colour = base_color), 
                        strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(colour = base_color), 
                        ...)
  th <- ggplot2::`%+replace%`(th1, th2)
  th
}





clusterData <- function (data, method = "hclust", rows = NULL, cols = NULL, 
          inf = NULL, na = NULL, ...) 
{

  if (is.null(rows)) 
    rows <- nrow(data) < 1000
  if (is.null(cols)) 
    cols <- ncol(data) < 1000
  if (is.null(na)) 
    na <- ""
  if (is.null(inf)) 
    inf <- ""
  if (any(!is.finite(nat(data)))) {
    if (naf(inf == FALSE)) 
      data <- subInf(data)
    if (is.na(inf)) 
      data[!is.finite(data)] <- inf
  }
  if (any(is.na(data))) {
    if (is.numeric(na)) 
      data[is.na(data)] <- na
  }
  tmp <- list(rows = data, cols = t(data))
  if (any(c(rows,cols))){
    res <- lapply(tmp[c(rows, cols)], function(data) {
      if (na == "omit") 
        data <- na.omit(data)
      clust <- NULL
      if (method == "hclust") 
        clust <- dendsort::dendsort(stats::hclust(stats::dist(data)))
      clust
    })
  } else {
    res <- NULL
  }

  res
}


ggseabar_bw <- function(results, x = NULL, y = NULL, label = NULL, sort_by = NULL, sort_abs = FALSE, top_by = NULL, top_n = 20, top_n_up = NULL, top_n_down = NULL, labsize = 12, reverse = FALSE, ...){
  
  # Input
  results <- as.data.frame(results)
  
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  label <- rlang::enquo(label)
  sort_by <- rlang::enquo(sort_by)
  top_by <- rlang::enquo(top_by)
  
  ordered_guess <- function(x, ref){
    guess <- grep(paste0(x, collapse = "|"), ref, value = TRUE)
    guess <- sort(factor(guess, ordered = TRUE, levels = x))
    guess[1]
  }
  
  if (rlang::quo_is_null(x)){
    x <- ordered_guess(c("padj","p.adjust","padjust","fdr", "FDR", "pvalue","pval"), colnames(results))
    x <- rlang::sym(as.character(x))
  }
  
  if (rlang::quo_is_null(y)){
    y <- ordered_guess( c("pathway", "Pathway", "description", "Description", "id", "ID", "geneset"), colnames(results))
    y <- rlang::sym(as.character(y))
  }
  
  if (rlang::quo_is_null(label)) label <- y
  if (rlang::quo_is_null(sort_by)) sort_by <- x
  if (rlang::quo_is_null(top_by)) top_by <- sort_by
  
  # Processing
  if (is.null(top_n_up)){
    results <- dplyr::arrange(results, dplyr::desc(abs(!!top_by)))
    results <- results[1:top_n,, drop = FALSE]
  } else {
    results <- dplyr::arrange(results, dplyr::desc(!!top_by))
    if (top_n_up > 0) ix.up <- 1:top_n_up else ix.up <- c()
    if (top_n_down > 0) ix.down <- (nrow(results) - top_n_down + 1):nrow(results) else ix.down <- c()
    ix <- unique(c(ix.up, ix.down))
    results <- results[ix,, drop = FALSE]
  }
  
  if (sort_abs == TRUE){
    results <- dplyr::arrange(results, dplyr::desc(abs(!!sort_by)))
  } else {
    results <- dplyr::arrange(results, dplyr::desc(!!sort_by))
  }
  
  if (reverse == TRUE){
    results[[rlang::as_name(y)]] <- factor(results[[rlang::as_name(y)]], ordered = TRUE, levels = results[[rlang::as_name(y)]][!duplicated(results[[rlang::as_name(y)]])])
  } else {
    results[[rlang::as_name(y)]] <- factor(results[[rlang::as_name(y)]], ordered = TRUE, levels = rev(results[[rlang::as_name(y)]][!duplicated(results[[rlang::as_name(y)]])]))
  }
  
  
  gg <- ggplot2::ggplot(results, aes(x = !!x, y = !!y, label = !!label, ...)) + ggplot2::geom_bar(stat = "identity", color = "black", fill = NA)
  gg <- gg + geom_vline(xintercept = 0, color = "grey70")
  
  res_left <- subset(results, results[[rlang::as_name(x)]] < 0)
  if (nrow(res_left) > 0) gg <- gg + ggplot2::geom_text(size = labsize/ggplot2::.pt, data = res_left, nudge_x = 0.25, hjust = 0)
  res_right <- subset(results, results[[rlang::as_name(x)]] > 0)
  if (nrow(res_right) > 0) gg <- gg + ggplot2::geom_text(size = labsize/ggplot2::.pt, data = res_right, nudge_x = -0.25, hjust = 1)
  
  lims <- ceiling(abs(range(results[[rlang::as_name(x)]])))
  gg <- gg + ylab("")
  gg <- gg + theme(axis.text.y = element_blank(),
                   panel.border = element_blank(),
                   plot.background = element_blank(),
                   axis.ticks.y = element_blank(),
                   plot.title = element_text(hjust = 0.5))
  
  
  gg <- gg + xlim(c(-max(abs(lims)), max(abs(lims))))
  
  gg
}




cxheatmap <- function (data, rowdf = NULL, coldf = NULL, scale = FALSE, cluster_rows = NULL, 
                       cluster_cols = NULL, rowdf_side = "left", coldf_side = "top", 
                       rowdf_legend = TRUE, coldf_legend = TRUE, legend_border = "black", 
                       anno_border = "black", fontsize = 12, rowcex = NULL, colcex = 1, 
                       rownames_width = 0.3, colnames_width = 0.3, heatpal = NULL, 
                       border = NULL, title = NULL, colors = NULL, inf = F, na = 0, 
                       mat = NULL, markoob = FALSE, markshape = 4, marksize = NULL, 
                       na_col = "grey", maxchar = 35, ...){
  
  datacall <- substitute(data)
  if (is.null(title)) {
    if (grepl("row|col", scale, ignore.case = TRUE)) {
      title <- "z-score"
    }
    else {
      title <- deparse1(datacall)
    }
  }
  if (title == FALSE) 
    title <- " "
  heatdata <- eval(datacall, envir = parent.frame())
  heatdata <- data.matrix(heatdata)
  heatdata <- matScale(heatdata, rows = grepl("row", scale, 
                                              ignore.case = TRUE), cols = grepl("col", scale, ignore.case = TRUE))
  clust <- clusterData(heatdata, rows = cluster_rows, cols = cluster_cols, 
                       inf = inf, na = na)
  if (is.null(clust$rows)) 
    clust$rows <- FALSE
  if (is.null(clust$cols)) 
    clust$cols <- FALSE
  if (is.null(heatpal)) {
    heatpal_colors <- getColorScale(heatdata)
    heatpal <- circlize::colorRamp2(breaks = heatpal_colors, 
                                    colors = names(heatpal_colors))
  }
  docol <- setdiff(unlist(lapply(list(coldf, rowdf), colnames)), 
                   names(colors))
  addcol <- NULL
  if (length(docol) > 0) {
    if (!is.null(coldf)) 
      all <- coldf
    if (!is.null(rowdf)) 
      all <- rowdf
    if (!is.null(coldf) & !is.null(rowdf)) 
      all <- dplyr::full_join(coldf, rowdf, by = character())
    addcol <- getColors(all[, docol, drop = FALSE])
  }
  colors <- c(colors, addcol)
  colors <- lapply(colors, function(tmp) tmp[!is.na(tmp) & 
                                               !is.na(names(tmp))])
  if (is.null(border)) {
    if (nrow(heatdata) < 100 & ncol(heatdata) < 100) {
      border <- grid::gpar(col = rgb(1, 1, 1), lwd = grid::unit(1, 
                                                                "pt"))
    }
    else {
      border <- grid::gpar(col = NA)
    }
  }
  else {
    if (any(border == TRUE)) {
      border <- grid::gpar(col = rgb(1, 1, 1), lwd = grid::unit(1, 
                                                                "pt"))
    }
    else if (length(border) > 1) {
      border <- grid::gpar(col = border[is.na(as.numeric(border))], 
                           lwd = grid::unit(as.numeric(border)[!is.na(as.numeric(border))], 
                                            "pt"))
    }
    else {
      border <- grid::gpar(col = NA)
    }
  }
  legend_params <- list(title_gp = grid::gpar(fontsize = fontsize, 
                                              fontface = "bold"), legend_height = grid::unit(0.2, "npc"), 
                        border = legend_border, labels_gp = grid::gpar(fontsize = fontsize))
  rowAnn <- NULL
  if (!is.null(rowdf)) 
    rowAnn <- getCXanno(df = rowdf[rownames(heatdata), , 
                                   drop = FALSE], colors = colors, anno_border = anno_border, 
                        side = rowdf_side, legend = rowdf_legend, legend_params = legend_params)
  colAnn <- NULL
  if (!is.null(coldf)) 
    colAnn <- getCXanno(coldf[colnames(heatdata), , drop = FALSE], 
                        colors = colors, anno_border = anno_border, side = coldf_side, 
                        legend = coldf_legend, legend_params = legend_params)
  if (markoob == TRUE & is.null(mat)) {
    mat <- matrix(data = FALSE, nrow = nrow(heatdata), ncol = ncol(heatdata), 
                  dimnames = dimnames(heatdata))
    mat[heatdata < min(heatpal_colors)] <- TRUE
    mat[heatdata > max(heatpal_colors)] <- TRUE
  }
  if (is.null(marksize)) 
    marksize <- fontsize * 0.6
  cellFUN <- NULL
  if (!is.null(mat)) {
    if (!is.logical(mat)) 
      stop("'Mat' must be a logical indicator of whether cells should be marked!")
    cellmat <- mat[rownames(heatdata), colnames(heatdata)]
    cellFUN <- function(j, i, x, y, width, height, fill) {
      if (naf(cellmat[i, j] == TRUE)) {
        grid::grid.points(x, y, pch = markshape, size = unit(marksize, 
                                                             "pt"))
      }
    }
  }
  dimnames(heatdata) <- lapply(dimnames(heatdata), function(x) cutstr(x, 
                                                                      maxchar = maxchar))
  if (is.null(rowcex) & nrow(heatdata) > 10 * ncol(heatdata)) 
    rowcex <- 1/log10(nrow(heatdata))
  if (is.null(rowcex)) 
    rowcex <- 1
  if (is.null(colcex)) 
    colcex <- 1
  hm <- ComplexHeatmap::Heatmap(name = title, matrix = heatdata, 
                                row_names_max_width = grid::unit(rownames_width, "npc"), 
                                column_names_max_height = grid::unit(colnames_width, 
                                                                     "npc"), column_title_gp = grid::gpar(fontsize = fontsize, 
                                                                                                          fontface = "bold"), rect_gp = border, na_col = na_col, 
                                left_annotation = rowAnn, top_annotation = colAnn, row_names_gp = grid::gpar(fontsize = fontsize * 
                                                                                                               rowcex), column_names_gp = grid::gpar(fontsize = fontsize * 
                                                                                                                                                       colcex), col = heatpal, heatmap_legend_param = legend_params, 
                                cluster_rows = clust$rows, cluster_columns = clust$cols, 
                                cell_fun = cellFUN, ...)
  hm
}





