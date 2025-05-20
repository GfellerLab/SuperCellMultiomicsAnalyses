CoveragePlot <- function (seurat, seurat.mc, metacell.col = NULL, sc.col = NULL, 
          dims = c(1, 2), reduction = "wnn.umap", mc.color = NULL, 
          sc.color = NULL, alpha = 1, pt_size = 0, metric = "size", 
          continuous_metric = F) 
{
  if (!is.null(seurat.mc@misc$membership)) {
    membership <- seurat.mc@misc$membership
  }
  else {
    if (!is.null(seurat.mc@misc$gamma)) {
      gamma <- seurat.mc@misc$gamma
    }
    else {
      gamma <- floor(ncol(seurat)/ncol(seurat.mc))
    }
    membership <- igraph::cut_at(seurat.mc@misc$metacells_hierarchy, 
                                 no = floor(ncol(seurat)/gamma))
  }
  seurat$Metacell <- membership
  seuratCoord <- Embeddings(seurat[[reduction]])
  seuratCoordMetacell <- cbind(seuratCoord, membership)
  centroids <- stats::aggregate(seuratCoord ~ membership, seuratCoord, 
                                mean)
  rownames(centroids) <- centroids[, 1]
  matches <- unlist(regmatches(colnames(seurat.mc), gregexpr("[[:digit:]]+", 
                                                             colnames(seurat.mc))))
  centroids <- centroids[matches, ]
  centroids[[metric]] <- seurat.mc[[metric]][, 1]
  if (is.null(metacell.col)) {
    metacell.col <- "SuperCell_MC"
    centroids[[metacell.col]] <- rep("red", length(seurat.mc$size))
  }
  else {
    centroids[[metacell.col]] <- seurat.mc[[metacell.col]][, 
                                                           1]
    print(head(centroids))
  }
  if (!is.null(sc.col)) {
    seuratCoord <- data.frame(seuratCoord)
    seuratCoord[[sc.col]] <- seurat[[sc.col]][, 1]
    p <- ggplot2::ggplot(seuratCoord, aes_string(colnames(seuratCoord)[dims[1]], 
                                                 colnames(seuratCoord)[dims[2]], color = sc.col)) + 
      ggplot2::stat_density_2d(alpha = alpha )
  }
  else {
    p <- ggplot2::ggplot(data.frame(seuratCoord), aes_string(colnames(seuratCoord)[dims[1]], 
                                                             colnames(seuratCoord)[dims[2]])) + ggplot2::stat_density_2d(alpha = alpha,
                                                                                                                    color = "grey")
  }
  if (!continuous_metric) {
    p <- p + ggplot2::geom_point(data = centroids, aes_string(colnames(centroids)[1 + 
                                                                                    dims[1]], colnames(centroids)[1 + dims[2]], fill = metacell.col, 
                                                              size = metric), colour = "black", pch = 21)
  }
  else {
    p <- p + ggplot2::geom_point(data = centroids, aes_string(colnames(centroids)[1 + 
                                                                                    dims[1]], colnames(centroids)[1 + dims[2]], fill = metric), 
                                 colour = "black", pch = 21, size = 2)
  }
  if (!is.null(metacell.col) & !is.null(sc.col) & !is.null(mc.color)) {
    if (metacell.col == sc.col & !is.null(mc.color)) {
      sc.color = mc.color
    }
  }
  if (!is.null(mc.color) & !continuous_metric) {
    p <- p + ggplot2::scale_fill_manual(values = mc.color)
  }
  if (!is.null(sc.color)) {
    p <- p + ggplot2::scale_color_manual(values = sc.color)
  }
  p <- p + ggplot2::theme_classic() + guides(fill = guide_legend(override.aes = list(size = 3.5)))
  return(p)
}