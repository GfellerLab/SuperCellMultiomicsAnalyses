## compactness and separation python functions need to be loaded with reticulate
metacellCompactnessSeparation <- function(sc.seurat,SC_seurat,preprocessing_method = "pca",dims =NULL) {
  if (is.null(dims)) {
    dims <- c(1:dim(sc.seurat@reductions[[preprocessing_method]]@cell.embeddings)[2])
  }
  memberships_without_outliers <- na.exclude(SC_seurat@misc$membership) # remove MetaCell2 outliers
  components = sc.seurat@reductions[[preprocessing_method]]@cell.embeddings[names(memberships_without_outliers),dims] 
  membership_df = data.frame(single_cell_names = names(memberships_without_outliers), 
                             SuperCell = sapply(memberships_without_outliers,function(i) paste0("SC_",i)))
  compactness_values = compactness(components, membership_df, label='SuperCell')
  separation_values = separation(components, membership_df, label='SuperCell')
  SC_seurat[[paste0("compactness_",preprocessing_method)]] <- compactness_values[paste0("SC_",colnames(SC_seurat)),"compactness"]
  SC_seurat[[paste0("separation_",preprocessing_method)]] <- separation_values[paste0("SC_",colnames(SC_seurat)),"separation"]
  
  return(SC_seurat)
}

RunDiffusionMap <- function(seurat, reduction = "pca", dims = c(1:50), dmDims = 10) {
  lowDimensionEmbeddings = seurat@reductions[[reduction]]@cell.embeddings[,dims] 
  cells <- colnames(seurat)
  dm <- PalantirDiffusionMap(lowDimensionEmbeddings,dims = dmDims) 
  rownames(dm) <- colnames(seurat)
  colnames(dm) <- c(1:length(colnames(dm)))
  seurat[["diffmap"]] <- CreateDimReducObject(embeddings = as.matrix(dm),
                                              assay = seurat@reductions[[reduction]]@assay.used,
                                              key = "DM_")
  return(seurat)
}


metacellSilhouetteSeurat <- function(seurat.mc,seurat.sc,reduction.name.sc,reduction,single.cells.res = F,dims=NULL) {
  if (is.null(dims)) {
    dims <- c(1:dim(seurat.sc@reductions[[reduction.name.sc]]@cell.embeddings)[2])
  }
  memberships_without_outliers <- na.exclude(seurat.mc@misc$membership) # remove MetaCell2 outliers
  
  dim.red = seurat.sc@reductions[[reduction.name.sc]]@cell.embeddings[names(memberships_without_outliers),dims]
  distances <- dist(dim.red)
  sil.res <- cluster::silhouette(memberships_without_outliers,distances)
  sil.res.summary <- summary(sil.res)
  seurat.mc@meta.data[names(sil.res.summary$clus.avg.widths),paste0("ASW_",reduction.name.sc)] <- sil.res.summary$clus.avg.widths
  if (single.cells.res) {
    seurat.mc@misc[[paste0("ASW_",reduction.name.sc)]] <- sil.res[,3]
  }
  return(seurat.mc)
}

metacellCompactnessSeurat <- function(seurat.mc,seurat.sc,reduction.name.sc,dims=NULL) {
  if (is.null(dims)) {
    dims <- c(1:dim(seurat.sc@reductions[[reduction.name.sc]]@cell.embeddings)[2])
  }
  
  memberships_without_outliers <- na.exclude(seurat.mc@misc$membership) # remove MetaCell2 outliers
  
  dim.red = seurat.sc@reductions[[reduction.name.sc]]@cell.embeddings[names(memberships_without_outliers),dims]
  centroids <- aggregate(dim.red,by = list(metacell =memberships_without_outliers),FUN = var)
  compactness <- apply(centroids[,-1], 1, mean)
  seurat.mc@meta.data[,paste0("compactness_",reduction.name.sc)] <- compactness
  return(seurat.mc)
}



mc_ASW <- function(seurat.mc,
                   seurat.sc,
                   reduction.name.sc,
                   reduction,
                   cell.membership,
                   group.label = "membership",
                   single.cells.res = F,dims=NULL) {
  if (is.null(dims)) {
    dims <- c(1:dim(seurat.sc@reductions[[reduction.name.sc]]@cell.embeddings)[2])
  }
  memberships_without_outliers <- na.exclude(cell.membership)
  membership_vector <- memberships_without_outliers[, group.label]
  names(membership_vector) <- rownames(memberships_without_outliers)
  #sc.reduction = sc.reduction[names(membership_vector), dims]
  
  dim.red = seurat.sc@reductions[[reduction.name.sc]]@cell.embeddings[names(membership_vector),dims]
  distances <- dist(dim.red)
  print(head(membership_vector))
  
  membership_vector_int <- as.integer(stringr::str_split_fixed(string = membership_vector,
                                                               pattern = "_",n = 2)[,2])
  print(head(membership_vector_int))
  prefix <- unique(stringr::str_split_fixed(membership_vector,
                                     "_",2)[,1])
  sil.res <- cluster::silhouette(membership_vector_int,distances)
  sil.res.summary <- summary(sil.res)
  metacell.asw <- sil.res.summary$clus.avg.widths
  names(metacell.asw) <- paste0(prefix,"_",names(metacell.asw))
  return(metacell.asw)
}