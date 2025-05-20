#' Weighted mean function
#' @export
#'
weighted.mean.fxn <- function (x, weights) {
  rowSums(x %*% Matrix::Diagonal(x = weights))/sum(weights)
}

#' Multimodal markers for all identity classes
#'
#' Function identifying markers in 2 modalities using FindAllMarkers.SuperCell
#'
#'
#' @param object Seurat metacell object with a metadata size column
#' @param assay1 Seurat assay to consider as the first modality
#' @param assay2 Seurat assay to consider as the second modality
#' @param test.use (default is "survey_weighted_t") Test to use in FindAllMarkers.SuperCell
#' @param fc.name1 fc.name parameter that FindAllMarkers.SuperCell will consider for the first modality
#' @param fc.name2 fc.name parameter that FindAllMarkers.SuperCell will consider for the second modality
#' @param logfc.threshold1 logfc.threshold parameter that FindAllMarkers.SuperCell will consider for the first modality
#' @param logfc.threshold2 logfc.threshold parameter that FindAllMarkers.SuperCell will consider for the second modality
#' @param mean.fxn1 mean.fxn parameter that FindAllMarkers.SuperCell will consider for the first modality
#' @param mean.fxn2 mean.fxn parameter that FindAllMarkers.SuperCell will consider for the second modality
#' @inheritParams SuperCellMultiomics::FindAllMarkers.SuperCell
#' @return A list of 3 data.frames: (i) a data frame containing a ranked list of putative multimodal markers as rows, and associated statistics as columns (p-values, t statistic, degree of freedom (df) of the weighted t-test).
#' avg_logFC, pct.1, pc.2 are computed taking into account metacell sizes, (ii) a data frame containing the output of FindAllMarkers.SuperCell run on the first modality, and (iii) a data frame containing the output of FindAllMarkers.SuperCell run onthe second modality.
#' @import Seurat
#' @export

FindMultimodalMarkers.SuperCell <- function(seurat.obj, group.by = 'celltype',
                                            assay1 = 'RNA', assay2 = "chromvar",
                                            min.cells.feature = 0, min.cells.group = 0, min.pct1 = 0.2,min.pct2 = 0,
                                            padj.cutoff = 0.05, base = 2, only.pos = T, return.thresh = 1,
                                            logfc.threshold1 = 0.1, logfc.threshold2 = 0.1,
                                            test.use = "survey_weighted_t",
                                            fc.name1 = "avg_log2FC", fc.name2 = "avg_diff",
                                            conversion = NULL,
                                            keep.only.best = FALSE,
                                            mean.fxn1 = NULL, mean.fxn2 = weighted.mean.fxn, ...) {
  
  features1 <-  rownames(seurat.obj[[assay1]])[rownames(seurat.obj[[assay1]])%in% conversion]
  print(head(features1))
  features2 <- rownames(seurat.obj[[assay2]])[rownames(seurat.obj[[assay2]])%in% names(conversion)]
  print(head(features2))

  n.features <- length(features1)
  print(n.features)
  
  Idents(seurat.obj) <- group.by
  markers_mod1 <- FindAllMarkers.SuperCell(
    object = seurat.obj, assay = assay1,
    min.cells.feature = min.cells.feature,
    min.cells.group = min.cells.group,
    base = base,
    features = features1,
    min.pct = min.pct1, fc.name = fc.name1,
    logfc.threshold = logfc.threshold1,
    only.pos = only.pos, mean.fxn = mean.fxn1, test.use = test.use,
    return.thresh = return.thresh, ...
  )
  
  markers_mod1$p_val_adj <- markers_mod1$p_val * n.features # correct bonferonni correction
  markers_mod1$p_val_adj[markers_mod1$p_val_adj >1] <- 1
  if (keep.only.best) {
    markers_mod1 <- data.frame(markers_mod1 %>%
      group_by(gene) %>%
      slice_min(p_val_adj,n = 1))
  }
 
  colnames(markers_mod1) <- paste0(assay1, ".", colnames(markers_mod1))
  
  markers_mod1$feature <- markers_mod1[, paste0(assay1, ".gene")]

  DefaultAssay(seurat.obj) <- assay2
  
  markers_mod2 <- FindAllMarkers.SuperCell(
    object = seurat.obj,
    assay = assay2,
    min.cells.feature = min.cells.feature,
    min.cells.group = min.cells.group,
    base = base,
    features = features2,
    min.pct = min.pct2,
    logfc.threshold = logfc.threshold2, test.use = test.use,
    only.pos = only.pos, fc.name = fc.name2, mean.fxn = weighted.mean.fxn, #SuperCellMultiomics:::weighted.mean.fxn
    return.thresh = return.thresh, ...
  )
  
  markers_mod2$p_val_adj <- markers_mod2$p_val * n.features # correct bonferonni correction
  markers_mod2$p_val_adj[markers_mod2$p_val_adj >1] <- 1
  
  
  if (keep.only.best) {
    markers_mod2 <- data.frame(markers_mod2 %>%
      group_by(gene) %>%
      slice_min(p_val_adj,n = 1))
  }
  
  colnames(markers_mod2) <- paste0(assay2, ".", colnames(markers_mod2))
  
  if (!is.null(conversion)) {
    markers_mod2$feature <- conversion[markers_mod2[, paste0(assay2, ".gene")]]
  }
  
  markers.all <- vector()
  for(celltype in unique(seurat.obj@meta.data[,group.by])){
    # print(celltype)
    ctmarkers_mod1 <- dplyr::filter(
      markers_mod1,
      !!dplyr::sym(paste0(assay1, ".cluster")) == celltype,
      !!dplyr::sym(paste0(assay1, ".p_val_adj")) < padj.cutoff,
      !!dplyr::sym(paste0(assay1, ".", fc.name1)) > 0) %>%
      dplyr::arrange(-!!dplyr::sym(paste0(assay1, ".", fc.name1)))
    
    ctmarkers_mod2 <- dplyr::filter(
      markers_mod2,
      !!dplyr::sym(paste0(assay2, ".cluster")) == celltype,
      !!dplyr::sym(paste0(assay2, ".p_val_adj"))  < padj.cutoff,
      !!dplyr::sym(paste0(assay2, ".", fc.name2)) > 0) %>%
      dplyr::arrange(-!!dplyr::sym(paste0(assay2, ".", fc.name2)))
    
    multimodal.markers <- dplyr::inner_join(
      x = ctmarkers_mod1,
      y = ctmarkers_mod2,
      by = "feature"
    )
    
    X <- as(GetAssayData(seurat.obj, assay = assay1), "dgCMatrix")
    y <- factor(ifelse(seurat.obj@meta.data[, group.by] == celltype, celltype, "other"))
    group.size <- as.numeric(table(y))
    n1n2 <- group.size * (ncol(X) - group.size)
    rank_res <- presto::rank_matrix(Matrix::t(X[multimodal.markers[, paste0(assay1, ".gene")],]))
    ustat <- presto:::compute_ustat(rank_res$X_ranked, y, n1n2, group.size)
    auc.rna <- t(ustat/n1n2)[,1]
    
    X <- as(GetAssayData(seurat.obj, assay = assay2), "dgCMatrix")
    rank_res <- presto::rank_matrix(Matrix::t(X[multimodal.markers[, paste0(assay2, ".gene")],]))
    ustat <- presto:::compute_ustat(rank_res$X_ranked, y, n1n2, group.size)
    auc.motif <- t(ustat/n1n2)[,1]
    
    multimodal.markers[, paste0(assay1, ".auc")] <- auc.rna
    multimodal.markers[, paste0(assay2, ".auc")]<- auc.motif #unlist(lapply(auc.res, function(x) x[,"motif.auc"]))
    multimodal.markers$auc.mean <- rowMeans(data.frame(auc.rna, auc.motif))
    
    if(test.use == "survey_weighted_t") {
    
    multimodal.markers[, paste0(assay1, ".rpb")] <- multimodal.markers[, paste0(assay1, ".t_value")] / sqrt(multimodal.markers[, paste0(assay1, ".t_value")]^2 + multimodal.markers[, paste0(assay1, ".df")])
    multimodal.markers[, paste0(assay2, ".rpb")] <- multimodal.markers[, paste0(assay2, ".t_value")] / sqrt(multimodal.markers[, paste0(assay2, ".t_value")]^2 + multimodal.markers[, paste0(assay2, ".df")])
    multimodal.markers$mean.rpb <- rowMeans(multimodal.markers[, c(paste0(assay1, ".rpb"), paste0(assay2, ".rpb"))])
    
    multimodal.markers <- dplyr::arrange(multimodal.markers, -mean.rpb)
    
    } else {
      multimodal.markers <- dplyr::arrange(multimodal.markers, auc.mean)
    }
    
    if(nrow(multimodal.markers) != 0){
      markers.all <- rbind(markers.all, multimodal.markers)
    }
  }
  
  return(
    setNames(
      list(
        markers.all,
        markers_mod1,
        markers_mod2
      ), nm = c("MultimodalMarkers", paste0("markers.", assay1), paste0("markers.", assay2))
    )
  )
}