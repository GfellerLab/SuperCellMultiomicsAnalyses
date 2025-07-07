
# Libraries ---------------------------------------------------------------
library(Seurat)
library(Signac)
library(getopt)
library(future)
library(presto)
library(BiocParallel)
library(SuperCellMultiomics)

# Increase max size limit
# options(future.globals.maxSize = 2000 * 1024^2)  # 2 GiB

# Parameters --------------------------------------------------------------

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds) or dataset name (from SeuratData package)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "scSeurat", "s", 1, "character", "Path to singlecell data"
), byrow=TRUE, ncol=5)

opt = getopt(spec)

dir.create(opt$outdir, recursive = T, showWarnings = F)
print(opt)

# Load single-cell data  ------------------------------
seurat.sc <- readRDS(opt$scSeurat)
DefaultAssay(seurat.sc) <- "RNA"
seurat.sc <- NormalizeData(seurat.sc)
seurat.sc <- FindVariableFeatures(seurat.sc, nfeatures = 2000)
var.genes <- VariableFeatures(seurat.sc)
saveRDS(var.genes, file = paste0(dirname(opt$scSeurat), "/hvgs.rds"))

# Load seurat object with metacells ---------------------------------------
if(opt$scSeurat != opt$inputSeurat){
  seurat.obj <- readRDS(opt$inputSeurat)
  use.weights <- T
}else{
  seurat.obj <- seurat.sc
  use.weights <- F
}

celltypes <- unique(seurat.obj$celltype)

############################################################################
############################################################################
###                                                                      ###
###                       CORRELATIONS COMPUTATION                       ###
###                                                                      ###
############################################################################
############################################################################

# Compute correlations between TF activities and gene expression ------------

# Extract TFs in common in the chromVAR and RNA assays
DefaultAssay(seurat.obj) <- "ATAC"
detectedMotifs <- rownames(seurat.obj[['chromvar']]) #motif
relatedTFs <- ConvertMotifID(seurat.obj, id = detectedMotifs) #tf

detectedMotifs <- detectedMotifs[which(relatedTFs %in% rownames(seurat.obj[["RNA"]]))]
relatedTFs <- relatedTFs[which(relatedTFs %in% rownames(seurat.obj[["RNA"]]))]

DefaultAssay(seurat.obj) <- 'RNA'
seurat.obj <- NormalizeData(seurat.obj)

chromVarCorrTable <- FeatureFeaturePlot.SuperCell(seurat.obj = seurat.obj,
                                                  is.normalized = T,
                                                  cluster = "celltype",
                                                  assays = c("RNA", "chromvar"),
                                                  feature.x = relatedTFs,
                                                  feature.y = detectedMotifs,
                                                  use.size = use.weights,
                                                  plot = F)

# Compute correlations between gene activities and gene expression ------------

# Extract genes in common between gene activity assays and RNA assay
GA.genes <- rownames(seurat.obj[["GA_ArchR"]])[rownames(seurat.obj[["GA_ArchR"]]) %in% rownames(seurat.obj[["RNA"]])]
common.genes <- var.genes[var.genes %in% GA.genes]

GA.ArchR.CorrTable <- FeatureFeaturePlot.SuperCell(seurat.obj = seurat.obj,
                                                   is.normalized = T,
                                                   cluster = "celltype",
                                                   assays = c("RNA", "GA_ArchR"),
                                                   feature.x = GA.genes,
                                                   feature.y = GA.genes,
                                                   use.size = use.weights,
                                                   plot = F)
ArchRGA.perCell <- cor(as.matrix(GetAssayData(seurat.obj, assay = "GA_ArchR", slot = "data")[GA.genes,]),
                       as.matrix(GetAssayData(seurat.obj, assay = "RNA", slot = "data")[GA.genes,]))


GA.genes <- rownames(seurat.obj[["GA_Signac"]])[rownames(seurat.obj[["GA_Signac"]]) %in% rownames(seurat.obj[["RNA"]])]
common.genes <- var.genes[var.genes %in% GA.genes]

GA.Signac.CorrTable <- FeatureFeaturePlot.SuperCell(seurat.obj = seurat.obj,
                                                    is.normalized = T,
                                                    cluster = "celltype",
                                                    assays = c("RNA", "GA_Signac"),
                                                    feature.x = GA.genes,
                                                    feature.y = GA.genes,
                                                    use.size = use.weights,
                                                    plot = F)
SignacGA.perCell <- cor(as.matrix(GetAssayData(seurat.obj, assay = "GA_Signac", slot = "data")[GA.genes,]),
                        as.matrix(GetAssayData(seurat.obj, assay = "RNA", slot = "data")[GA.genes,]))

cor.matrices <- list(chromVar.RNA = chromVarCorrTable,
                     ArchRGA.RNA = GA.ArchR.CorrTable,
                     SignacGA.RNA = GA.Signac.CorrTable,
                     ArchRGA.perCell = ArchRGA.perCell,
                     SignacGA.perCell = SignacGA.perCell)

saveRDS(cor.matrices, paste0(opt$outdir,"/CorrTables.rds"))

############################################################################
############################################################################
###                                                                      ###
###                          MULTIMODAL MARKERS                          ###
###                                                                      ###
############################################################################
############################################################################

DefaultAssay(seurat.obj) <- "ATAC"
if(!use.weights){
  # Extract top TFs at the single-cell level --------------------------------
  seurat.obj$size <- 1
  multimodalMarkersMethod = "weighted_t"
  markers.summary <- FindMultimodalMarkers.SuperCell(seurat.obj = seurat.obj, group.by = 'celltype',
                                                     assay1 = "RNA", assay2 = "chromvar",
                                                     min.cells.feature = 0, min.cells.group = 0,
                                                     min.pct = 0.01, padj.cutoff = 0.05, base = 2,
                                                     test.use = multimodalMarkersMethod,
                                                     only.pos = T, fc.name1 = "avg_log2FC",
                                                     fc.name2 = "avg_diff")

  saveRDS(markers.summary, file = paste0(opt$outdir, "/multimodalMarkers_ttest.rds"))


  # Get TopTFs using Seurat approach ----------------------------------------

  wilcox.rna <- presto:::wilcoxauc.Seurat(
    X = seurat.obj,
    group_by = 'celltype',
    assay = 'data',
    seurat_assay = 'RNA'
  )
  wilcox.motifs <- presto:::wilcoxauc.Seurat(
    X = seurat.obj,
    group_by = 'celltype',
    assay = 'data',
    seurat_assay = 'chromvar'
  )

  motif.names <- wilcox.motifs$feature
  colnames(wilcox.rna) <- paste0("RNA.", colnames(wilcox.rna))
  colnames(wilcox.motifs) <- paste0("chromvar.", colnames(wilcox.motifs))
  wilcox.rna$gene <- wilcox.rna$RNA.feature
  DefaultAssay(seurat.obj) <- "ATAC"
  wilcox.motifs$gene <- ConvertMotifID(seurat.obj, id = motif.names)

  saveRDS(wilcox.rna, file = paste0(opt$outdir, "wilcox_rna.rds"))
  saveRDS(wilcox.motifs, file = paste0(opt$outdir, "wilcox_motifs.rds"))

  topTFs.original <- function(markers_rna, markers_motifs, celltype, padj.cutoff = 0.05) {
    ctmarkers_rna <- dplyr::filter(
      markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0.1) %>%
      dplyr::arrange(-RNA.auc)
    ctmarkers_motif <- dplyr::filter(
      markers_motifs, chromvar.group == celltype, chromvar.padj < padj.cutoff, chromvar.logFC > 0.1) %>%
      dplyr::arrange(-chromvar.auc)
    top_tfs <- dplyr::inner_join(
      x = ctmarkers_rna[, c(2, 11, 6, 8)],
      y = ctmarkers_motif[, c(2, 1, 11, 6, 8)], by = "gene"
    )
    top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$chromvar.auc) / 2
    top_tfs <- dplyr::arrange(top_tfs, -avg_auc)
    return(top_tfs)
  }

  topTFs.wilcoxonAUC <- do.call(rbind, lapply(unique(seurat.obj$celltype), function(i) topTFs.original(markers_rna = wilcox.rna, markers_motifs = wilcox.motifs, celltype = i)))
  saveRDS(topTFs.wilcoxonAUC, paste0(opt$outdir,"/multimodalMarkers_wilcoxonAUC.rds"))


}else{
  for(multimodalMarkersMethod in c("nonWeigthed", "survey_weighted_t", "weighted_t")){
    if(multimodalMarkersMethod == "nonWeigthed"){
      seurat.obj$size <- 1
      test.method = "weighted_t"
    }else{
      test.method = multimodalMarkersMethod
    }

    markers.summary <- FindMultimodalMarkers.SuperCell(seurat.obj = seurat.obj, group.by = 'celltype',
                                                       assay1 = "RNA", assay2 = "chromvar",
                                                       min.cells.feature = 0, min.cells.group = 0,
                                                       min.pct = 0.01, padj.cutoff = 0.05, base = 2,
                                                       test.use = test.method,
                                                       only.pos = T,
                                                       fc.name1 = "avg_log2FC", fc.name2 = "avg_diff")

    saveRDS(markers.summary, paste0(opt$outdir,"/multimodalMarkers_", gsub("_", "",multimodalMarkersMethod), ".rds"))

  }


}
