FeatureFeaturePlot.SuperCell <- function (seurat.obj, feature.x, feature.y, cluster = NULL,
                                          assays = c("RNA", "RNA"), #nb.cores = 2,
                                          norm.margin.1 = 1, normalization.method.1 = "LogNormalize",
                                          norm.margin.2 = 1, normalization.method.2 = "LogNormalize",
                                          is.normalized = F, plot = F, color.use = NULL, use.size = T, add.key = T,...)
{
  Seurat::DefaultAssay(seurat.obj) <- assays[1]
  if (!is.normalized) {
    seurat.obj <- Seurat::NormalizeData(seurat.obj, normalization.method = normalization.method.1,
                                        margin = norm.margin.1)
  }
  fe1 <- Seurat::GetAssayData(seurat.obj, slot = "data", assay = assays[1])[feature.x, , drop= F ]
  feature.x <- paste0(gsub("_", "", tolower(assays[1])), "_", feature.x)
  rownames(fe1) <- feature.x
  
  Seurat::DefaultAssay(seurat.obj) <- assays[2]
  if (!is.normalized) {
    seurat.obj <- Seurat::NormalizeData(seurat.obj, normalization.method = normalization.method.2,
                                        margin = norm.margin.2)
  }
  fe2 <- Seurat::GetAssayData(seurat.obj, slot = "data", assay = assays[2])[feature.y, , drop= F]
  feature.y <- paste0(gsub("_", "", tolower(assays[2])), "_", feature.y)
  rownames(fe2) <- feature.y
  fe <- rbind(fe1, fe2)
  
  if (use.size) {
    sizes <- as.numeric(seurat.obj$size)
  } else {
    sizes <- rep(1, ncol(seurat.obj))
  }
  
  res <- list()
  cor.res <- lapply(1:length(feature.x), function(i) {
    tryCatch({
      weights::wtd.cor(fe[feature.x[i], ], fe[feature.y[i], ], weight = sizes)
    }, error = function(e) {
      data.frame(correlation = NA, p.value = NA)
    })
  })
  
  res[["w.cor"]] <- lapply(cor.res, function(x) x[,"correlation"])
  res[["w.pval"]] <- lapply(cor.res, function(x) x[, "p.value"])
  if(plot){
    seurat.obj$size <- sizes
    res[["p"]] <- lapply(1:length(feature.x), function(i){
      if(add.key){
        FeatureScatter.SuperCell(object = seurat.obj,
                                 size.by = "size",
                                 feature1 = feature.x[i],
                                 feature2 = feature.y[i],
                                 group.by = cluster,
                                 cols = color.use,
                                 plot.cor = plot,...)
      }else{
        FeatureScatter.SuperCell(object = seurat.obj,
                                 size.by = "size",
                                 feature1 = feature.x[i],
                                 feature2 = feature.y[i],
                                 group.by = cluster,
                                 cols = color.use,
                                 plot.cor = plot,...) +
          xlab(gsub(paste0(gsub("_", "", tolower(assays[1])),"_"), "", feature.x[i])) +
          ylab(gsub(paste0(gsub("_", "", tolower(assays[2])),"_"), "", feature.y[i]))
      }
    } )
    names(res[["p"]]) <- feature.x
  }
  
  w.cor <- data.frame(row.names = paste(feature.x, feature.y, sep = "_"),
                      feature1 = feature.x,
                      feature2 = feature.y,
                      w.cor = as.numeric(res$w.cor),
                      w.pval = as.numeric(res$w.pval),
                      w.qval = p.adjust(as.numeric(res$w.pval), method = "bonferroni"))
  
  if (plot) {
    return(list(cor.res = w.cor, plots = res$p))
  }else{
    return(w.cor)
  }
}