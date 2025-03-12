subsampleCellType <- function(seurat,
                              cellType,
                              n=20,
                              nfeatures = 2000,
                              cellTypeCol ="cell_line",
                              RNAnormalization = "LogNormalize",
                              seed = 1) {
  Idents(seurat) <- cellTypeCol
  subIdx = which(seurat[[cellTypeCol]][,1] == cellType)
  print(subIdx)
  popSize <- length(subIdx)
  set.seed(seed)
  subCells <- sample(subIdx,size = popSize- n)
  seurat <- subset(seurat,cells = subCells,invert = T)
  print(table(seurat[[cellTypeCol]]))
  
  #seurat <- NormalizeData(seurat)
  
  if (RNAnormalization == "SCT") {
    rnaAssay = "SCT"
    DefaultAssay(seurat) <- "RNA"
    seurat <- SCTransform(seurat, verbose = FALSE) 
    seurat <- RunPCA(seurat)
  } else {
    rnaAssay = "RNA"
    seurat <- NormalizeData(seurat, verbose = FALSE) %>% FindVariableFeatures(,nFeature = nfeatures) %>% ScaleData() %>% RunPCA()
  }
  
  
  if ("ADT" %in% names(seurat@assays)) {
    DefaultAssay(seurat) <- 'ADT'
    # we will use all ADT features for dimensional reduction
    # we set a dimensional reduction name to avoid overwriting the 
    VariableFeatures(seurat) <- rownames(seurat[["ADT"]])
    seurat <- NormalizeData(seurat, normalization.method = 'CLR', margin = 2) %>% 
      ScaleData() %>% RunPCA(reduction.name = 'apca')
  }
  
  if ("ATAC" %in% names(seurat@assays)) {
    DefaultAssay(seurat) <- 'ATAC'
    seurat <- RunTFIDF(seurat)
    seurat <- FindTopFeatures(seurat, min.cutoff = "q0")
    seurat <- RunSVD(seurat)
  }
  return(seurat)
}


findRarePop <- function(seuratWithRarePop,
                        cellTypeCol,
                        cellType,
                        kernel = T, 
                        k.knn = 30, 
                        assay = c("RNA"),
                        reduction = list("pca"),
                        graph.name = "nn",
                        dims = c(1:30)) {
  N.SC <- length(unique(seuratWithRarePop[[cellTypeCol]][,1])) -1
  seuratWithRarePop.mc.0 <- SCimplify_for_Seurat(seurat =seuratWithRarePop,
                                                 gamma = ncol(seuratWithRarePop)/N.SC,
                                                 graph.name = graph.name,
                                                 assay = assay,
                                                 reduction = reduction,
                                                 kernel = kernel, 
                                                 k.knn = k.knn, 
                                                 dims = dims,
                                                 return.seurat = T)
  notDetected <- T
  cluster <- seuratWithRarePop[[cellTypeCol]][,1]
  while (notDetected) {
    N.SC <- N.SC+1
    #print(N.SC)
    newGamma <- ncol(seuratWithRarePop)/N.SC
    #print(newGamma)
    SC <- SCimplify_for_Seurat(seurat =seuratWithRarePop,
                               seurat.mc = seuratWithRarePop.mc.0,
                               gamma = newGamma,
                               kernel = kernel, 
                               k.knn = k.knn, 
                               dims = dims,
                               return.seurat = F)
    
    SC[[cellTypeCol]] <- supercell_assign(clusters = cluster,supercell_membership = SC$membership,method = "absolute")
    SC[["purity"]] <- supercell_purity(clusters = cluster,supercell_membership = SC$membership)
    
    #SC$supercell_size <- as.numeric(table(SC$membership))
    notDetected <- ! cellType %in% unique(SC[[cellTypeCol]])
    
  }
  res <- list(k = N.SC ,
              gamma = length(cluster)/N.SC,
              prop = (SC$supercell_size[SC[[cellTypeCol]] == cellType]*SC$purity[SC[[cellTypeCol]] == cellType])/length(which(cluster == cellType)),
              purity = SC$purity[SC[[cellTypeCol]] == cellType] )
  
  return(res)
}