## compactness and separation python functions need to be loaded with reticulate
metacellCompactnessSeparation <- function(sc.seurat,SC_seurat,preprocessing_method = "pca",dims =c(1:50)) {
  components = sc.seurat@reductions[[preprocessing_method]]@cell.embeddings[,dims] 
  membership_df = data.frame(single_cell_names = names(SC_seurat@misc[["membership"]]), 
                             SuperCell = sapply(SC_seurat@misc[["membership"]],function(i) paste0("SC_",i)))
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



