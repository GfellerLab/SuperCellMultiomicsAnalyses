ProjectMetacellsUnimodalUMAP <- function(sobj.sc,
                                         sobj.mc,
                                         umap.sc,
                                         dim.red.sc,
                                         dims = NULL) {
  
  assay.used <- sobj.sc[[dim.red.sc]]@assay.used
  
  if (is.null(dims)) {
    dims <- 1:sobj.sc[[dim.red.sc]]@misc$model$norig_col #Be carefull won't be valid when 2:N has been used (eg lsi atac)
  }
  
  sobj.sc$Metacell <- sobj.mc@misc$membership
  
  
  ##Get metacell centroid in dim.red.sc
  dim.red.sc.coord <- Embeddings(bm[[dim.red.sc]])
  # dim.red.sc.coord <- cbind(dim.red.sc.coord,  
  #                           data.frame("membership" = sobj.mc@misc$membership))
  centroids <- stats::aggregate(dim.red.sc.coord ~  sobj.mc@misc$membership, 
                                dim.red.sc.coord, 
                                mean)
  centroids <-  as.matrix(centroids)[,-1]
  centroids <- centroids[,dims]
  rownames(centroids) <- colnames(sobj.mc)
  sobj.mc[[dim.red.sc]] <- CreateDimReducObject(embeddings =centroids,
                                                assay = assay.used,
                                                key = paste("mean_sc_",dim.red.sc))
  
  sobj.mc <- ProjectUMAP(
    query = sobj.mc,
    query.reduction = dim.red.sc,
    reference = bm,
    reduction.name = paste0("mean_",dim.red.sc,"_sc_proj_",umap.sc),
    reference.reduction = dim.red.sc,
    reduction.model = umap.sc
  )
  
  return(sobj.mc)
  
}



CoveragePlotMetacells <- function(sobj.sc,
                                  reduction.sc,
                                  sobj.mc,
                                  reduction.mc,
                                  group.by = NULL,
                                  cols =NULL,
                                  dims = c(1, 2)) {
  
  dimplot.sc <- DimPlot(sobj.sc,
                        dims = dims,
                        reduction = reduction.sc,
                        group.by = group.by,
                        cols = cols
  ) 
  
  embedding.keys.mc <- paste0(sobj.mc[[reduction.mc]]@key,dims)
  data.mc <- FetchData(sobj.mc,vars = c(embedding.keys.mc,"size",group.by))
  
  covplot <- dimplot.sc + ggplot2::geom_point(data = data.mc,
                                              aes_string(embedding.keys.mc[1],
                                                         embedding.keys.mc[2],
                                                         fill = group.by, 
                                                         size = "size"), 
                                              colour = "black", pch = 21) + 
    scale_fill_manual(values = cols) + 
    guides(fill = "none")
  
  return(covplot)
}


GetJointSingleCellsMetacells <- function(sobj.sc,
                                         sobj.mc,
                                         dimred_1 = "pca",
                                         dimred_2= "apca") {
  
  assay_1 <- sobj.sc[[dimred_1]]@assay.used
  assay_2 <- sobj.sc[[dimred_2]]@assay.used
  
  
  sobj.sc <- DietSeurat(sobj.sc,dimreducs = c(dimred_1,dimred_2),layers = "counts",
                        features = c(rownames(sobj.sc[[assay_1]])[1:2],
                                     rownames(sobj.sc[[assay_2]])[1:2]
                        ))
  sobj.mc <- DietSeurat(sobj.mc,features = c(rownames(sobj.sc[[assay_1]])[1:2],
                                             rownames(sobj.sc[[assay_2]])[1:2]
  ))
  
  
  ## dimred 1
  seuratCoord <- Embeddings(sobj.sc[[dimred_1]])
  
  centroids <- stats::aggregate(seuratCoord ~  sobj.mc@misc$membership, seuratCoord, 
                                mean)
  
  centroids <-  as.matrix(centroids)[,-1]
  rownames(centroids) <- colnames(sobj.mc)
  
  sobj.mc[[dimred_1]] <- CreateDimReducObject(embeddings =centroids,assay = sobj.sc[[dimred_1]]@assay.used)
  
  ## dimred 2
  seuratCoord <- Embeddings(sobj.sc[[dimred_2]])
  centroids <- stats::aggregate(seuratCoord ~  sobj.mc@misc$membership, seuratCoord, 
                                mean)
  
  centroids <-  as.matrix(centroids)[,-1]
  rownames(centroids) <- colnames(sobj.mc)
  
  sobj.mc[[dimred_2]] <- CreateDimReducObject(embeddings =centroids,assay = sobj.sc[[dimred_2]]@assay.used)
  joint.sobj <- merge(sobj.sc,sobj.mc)
  
  joint.sobj[[dimred_1]] <- merge(sobj.sc[[dimred_1]],sobj.mc[[dimred_1]])
  joint.sobj[[dimred_2]] <- merge(sobj.sc[[dimred_2]],sobj.mc[[dimred_2]])
  
  return(joint.sobj)
}

SingleCellsMetacellsWnnUMAP <- function(sobj.sc,
                                        sobj.mc,
                                        dimred_1 = "pca",
                                        dimred_2= "apca",
                                        dims.list = list(1:30, 1:18),
                                        dims.umap = c(1,2),
                                        cols = NULL,
                                        group.by = NULL) {
  
  joint.sobj <- GetJointSingleCellsMetacells(sobj.sc = sobj.sc,
                                             sobj.mc = sobj.mc,
                                             dimred_1 = dimred_1,
                                             dimred_2 =  dimred_2)
  
  joint.sobj <- FindMultiModalNeighbors(
    joint.sobj,
    reduction.list = list(dimred_1, dimred_2),
    dims.list = dims.list
  )
  
  joint.sobj <- RunUMAP(
    joint.sobj,
    nn.name = "weighted.nn",
    reduction.name = "wnn.umap"
  )
  singlecells <- colnames(joint.sobj)[is.na(joint.sobj$size)]
  metacells <- colnames(joint.sobj)[!is.na(joint.sobj$size)]
  
  
  sobj.sc[["sc_mc_wnn_umap"]] <- CreateDimReducObject(
    embeddings = Embeddings(joint.sobj[["wnn.umap"]])[singlecells,]
  )
  
  sobj.mc[["sc_mc_wnn_umap"]] <- CreateDimReducObject(
    embeddings = Embeddings(joint.sobj[["wnn.umap"]])[metacells,]
  ) 
  
  covplot <-CoveragePlotMetacells(sobj.sc = sobj.sc,
                                  sobj.mc = sobj.mc,
                                  reduction.sc = "sc_mc_wnn_umap",
                                  reduction.mc = "sc_mc_wnn_umap",
                                  cols = cols,
                                  group.by = group.by,
                                  dims = dims.umap)
  
  return(covplot)
  
}


