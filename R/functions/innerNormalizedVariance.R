inv <- function(x, logScale = F){
  if (logScale) {
    res <- log(var(x)/mean(x))
  } else {
    res <- var(x)/mean(x)
  }
  return(res)
}

innerNormVar = function(gene,seurat,memberships) {
  x = seurat@assays$RNA@counts[gene,]
  res <- aggregate(x=as.vector(x),by=list(metacell = memberships),FUN = "inv")[,"x"]
  return(res)
}


computeInnerNormVar <- function(seurat, memberships, genesToTest = NULL) {
  
  if (is.null(genesToTest)) {
    DefaultAssay(seurat) <- "RNA"
    seurat <- Seurat::FindVariableFeatures(seurat)
    genesToTest <- Seurat::VariableFeatures(seurat,assay = "RNA")
  }
  resTable <- sapply(X = genesToTest,FUN = innerNormVar, seurat,memberships = memberships)
  res <- apply(resTable, MARGIN = 1, FUN = quantile,probs = 0.95, na.rm = T)
  return(res)
}