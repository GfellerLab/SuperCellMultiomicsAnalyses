
crPLot <- function (seurat.mc, feature_x, feature_y, method = c("pearson", 
                                                                "kendall", "spearman"), assays = c("ACTIVITY", "RNA"), cluster = "celltype.l2", 
                    is.normalized = T, plot = T, color.use = NULL, use.size = T) 
{
  method <- match.arg(arg = method)
  Seurat::DefaultAssay(seurat.mc) <- assays[1]
  if (assays[1] == "RNA") {
    if (!is.normalized) {
      seurat.mc <- Seurat::NormalizeData(seurat.mc, normalization.method = "LogNormalize", 
                                         margin = 1)
    }
  }
  else {
    if (!is.normalized) {
      seurat.mc <- Seurat::NormalizeData(seurat.mc, normalization.method = "CLR", 
                                         margin = 2)
    }
  }
  fe1 <- Seurat::GetAssayData(seurat.mc, slot = "data", assay = assays[1])[feature_x, 
  ]
  feature_x <- paste0(tolower(assays[1]), "_", feature_x)
  rownames(fe1) <- feature_x
  Seurat::DefaultAssay(seurat.mc) <- assays[2]
  if (assays[2] == "ADT") {
    if (!is.normalized) {
      seurat.mc <- Seurat::NormalizeData(seurat.mc, normalization.method = "LogNormalize", 
                                         margin = 1)
    }
  }
  else {
    if (!is.normalized) {
      seurat.mc <- Seurat::NormalizeData(seurat.mc, normalization.method = "CLR", 
                                         margin = 2)
    }
  }
  fe2 <- Seurat::GetAssayData(seurat.mc, slot = "data", assay = assays[2])[feature_y, 
  ]
  rownames(fe2) <- feature_y
  fe <- rbind(fe1, fe2)
  if (use.size) {
    sizes <- as.numeric(seurat.mc$size)
  }
  else {
    sizes <- rep(1, length(ncol(seurat.mc)))
  }
  res <- supercell_FeatureFeaturePlot(fe, feature_x = feature_x, 
                                      feature_y = feature_y, method = method, supercell_size = sizes, 
                                      cluster = seurat.mc[[cluster]][, 1], color.use = color.use, 
                                      combine = F)
  w.cor <- data.frame(features = names(res$w.cor), w.cor = as.numeric(res$w.cor))
  if (plot) {
    for (i in names(res$p)) {
      plot(res$p[[i]])
    }
  }
  return(list(w.cor,res))
}



supercell_FeatureFeaturePlot_single <- function(fe_x,
                                                fe_y,
                                                feature_x_name,
                                                feature_y_name,
                                                method = c("pearson", "kendall", "spearman"),
                                                supercell_size = NULL,
                                                clusters = NULL,
                                                color.use = NULL,
                                                x.max = NULL,
                                                y.max = NULL,
                                                pt.size = 1,
                                                x.min = NULL,
                                                y.min = NULL,
                                                alpha = 0.9){
  
  
  method <- match.arg(arg = method)
  N.SC <- length(fe_x)
  
  plot.df <- data.frame(x = fe_x,
                        y = fe_y,
                        identity = factor(clusters),
                        size = supercell_size)
  
  
  membership <- rep(1:N.SC, plot.df$size)
  
  crt               <- stats::cor.test(plot.df$x[membership], plot.df$y[membership],method = method)
  w.cor             <- unname(crt$estimate)
  w.pval            <- unname(crt$p.value)
  
  if(is.null(x.max)) x.max <- NA
  if(is.null(y.max)) y.max <- NA
  
  if(is.null(x.min)) x.min <- NA
  if(is.null(y.min)) y.min <- NA
  
  if(length(unique(plot.df$size)) == 1){
    
    g <- ggplot2::ggplot(plot.df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$identity)) +
      ggplot2::geom_point(alpha = alpha, size = pt.size)
  } else {
    g <- ggplot2::ggplot(plot.df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$identity, size = .data$size)) +
      ggplot2::geom_point(alpha = alpha)
  }
  
  g <- g +
    ggplot2::scale_x_continuous(limits = c(x.min, x.max)) +
    ggplot2::scale_y_continuous(limits = c(y.min, y.max)) +
    ggplot2::theme_classic() + ggplot2::theme(aspect.ratio = 1) + #, legend.position = "none"
    ggplot2::labs(x = feature_x_name,
                  y = feature_y_name,
                  title = paste0("w.cor = ", signif(w.cor, 2))) #remove pval in the title
  
  
  if(!is.null(color.use)){
    g <- g + ggplot2::scale_color_manual(values = color.use)
  }
  
  res <- list(g = g, w.cor = w.cor, w.pval = w.pval)
  return(res)
}


supercell_FeatureFeaturePlot <- function(fe,
                                         feature_x,
                                         feature_y,
                                         method = c("pearson", "kendall", "spearman"),
                                         supercell_size = NULL,
                                         clusters = NULL,
                                         color.use = NULL,
                                         idents = NULL,
                                         pt.size = 1,
                                         alpha = 0.9,
                                         x.max = NULL,
                                         y.max = NULL,
                                         same.x.lims = FALSE,
                                         same.y.lims = FALSE,
                                         ncol = NULL,
                                         combine = TRUE,
                                         sort.by.corr = TRUE){
  
  method <- match.arg(arg = method)
  N.SC <- ncol(fe) # number of super-cells
  
  if(is.null(clusters)) clusters <- 1
  
  if((length(clusters) != N.SC) & length(clusters) != 1){
    stop(paste0("clusters has to be a vector of the same lenght as fe1 (", N.SC, ") or 1, not ", length(clusters)))
  }
  if(length(clusters) == 1){
    clusters <- rep(clusters, N.SC)
  }
  
  
  if(is.null(supercell_size)) supercell_size <- rep(1, N.SC)
  
  if((length(supercell_size) != N.SC) & length(supercell_size) != 1){
    stop(paste0("supercell_size has to be a vector of the same lenght as fe1 (", N.SC, ") or 1, not ", length(supercell_size)))
  }
  if(length(supercell_size) == 1){
    supercell_size <- rep(supercell_size, N.SC)
  }
  
  
  if(is.null(idents)) idents <- sort(unique(clusters))
  
  ids.keep.idents <- which(clusters %in% idents)
  
  
  if(!is.null(color.use)){
    if(length(color.use) < length(idents)){
      warning(paste0("Length of color.use (", length(color.use), ") is smaller than number of idents (",
                     length(idents),"), color.use will not be used"))
      color.use <- NULL
    }
  }
  
  
  if(length(feature_x) != length(feature_y)){
    if(length(feature_x) == 1){
      feature_x <- rep(feature_x, length(feature_y))
    } else {
      if(length(feature_y) == 1){
        feature_y <- rep(feature_y, length(feature_x))
      } else{
        stop("Vectors feature_x and feature_y need to have the same length or one of them has to be a vector of length 1")
      }
    }
  }
  
  # keep features that are present in the feature expression dataset
  features.set_x <- feature_x[feature_x %in% rownames(fe) & feature_y %in% rownames(fe)]
  features.set_y <- feature_y[feature_x %in% rownames(fe) & feature_y %in% rownames(fe)]
  
  
  if(same.y.lims & is.null(x = y.max)){
    y.max <- max(fe[features.set_y, ids.keep.idents])
  }
  
  if(same.x.lims & is.null(x.max)){
    x.max <- max(fe[features.set_x, ids.keep.idents])
  }
  
  p.list <- list()
  
  
  
  for(i in 1:length(features.set_x)){
    
    features.i <- paste(features.set_x[i], features.set_y[i], sep = "_")
    
    p.list[[features.i]] <- supercell_FeatureFeaturePlot_single(fe_x = fe[features.set_x[i], ids.keep.idents],
                                                                fe_y = fe[features.set_y[i], ids.keep.idents],
                                                                method = method,
                                                                feature_x_name = features.set_x[i],
                                                                feature_y_name = features.set_y[i],
                                                                supercell_size = supercell_size[ids.keep.idents],
                                                                clusters = clusters[ids.keep.idents],
                                                                color.use = color.use,
                                                                x.max = x.max,
                                                                y.max = y.max,
                                                                x.min =NULL,
                                                                y.min =NULL,
                                                                pt.size = pt.size,
                                                                alpha = alpha)
  }
  
  
  
  if(sort.by.corr){ # sort plots by absolute value of correkation
    p.list <- p.list[names(sort(abs(unlist(lapply(p.list, FUN = function(x){x$w.cor}))), decreasing = T, na.last = T))]
  }
  
  w.cor.list <- lapply(p.list, FUN = function(x){x$w.cor})
  w.pval.list <- lapply(p.list, FUN = function(x){x$w.pval})
  p.list <-lapply(p.list, FUN = function(x){x$g})
  
  if(combine) {
    p.list <- patchwork::wrap_plots(p.list, ncol = ncol, guides = "collect")
  }
  return(list(p = p.list, w.cor = w.cor.list, w.pval = w.pval.list))
}