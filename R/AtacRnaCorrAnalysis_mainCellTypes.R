
# Libraries ---------------------------------------------------------------
library(Seurat)
library(Signac)
library(getopt)
library(future)
library(presto)
library(BiocParallel)
library(SuperCell)

# Increase max size limit
# options(future.globals.maxSize = 2000 * 1024^2)  # 2 GiB
# Functions ---------------------------------------------------------------
source("R/functions/SuperCellMultiomics_functions2.R")


# Parameters --------------------------------------------------------------

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds) or dataset name (from SeuratData package)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "scSeurat", "s", 1, "character", "Path to singlecell data"
), byrow=TRUE, ncol=5)

opt = getopt(spec)

# machine_path <- "" # "/mnt/curnagl/"
# dataset <- "pbmcMultiome" #pbmcMultiome hspcMultiomePersad
# data_path <- paste0(machine_path, "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/SuperCellMultiomicsAnalyses/manuscript_version/output/", dataset, "/")
# opt <- list()
# # opt$inputSeurat <- paste0("/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/SuperCellMultiomicsAnalyses/manuscript_version/output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds")
# opt$inputSeurat <- paste0("/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/SuperCellMultiomicsAnalyses/manuscript_version/output/pbmcMultiome/SuperCellMulti/g20/seurat.multiome.activities.rds") #seurat.multiome.ArchRGA.mc.rds"
# opt$scSeurat <- paste0("/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/SuperCellMultiomicsAnalyses/manuscript_version/output/pbmcMultiome/singlecells_analysis/seurat.multiome.activities.rds")
# opt$outdir <- paste0(data_path, "/SuperCellMulti/g20")

# if(is.null(opt$nWorkers)){
#   opt$nWorkers = parallel::detectCores()-4
# }
# future::plan("multisession", workers = opt$nWorkers)

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

seurat.obj <- seurat.obj[, !seurat.obj$celltype %in% c("Plasma", "HSPC")]
# seurat.obj$celltype_main <- plyr::revalue(seurat.obj$celltype, c("B Interm" = "Bcells", "B Mem" = "Bcells", "B Naive" = "Bcells",
#                                                                "CD14 Mono" = "Monocytes", "CD16 Mono" = "Monocytes",
#                                                                "Treg"= "CD4_Tcells","CD4 Naive" = "CD4_Tcells", "CD4 Mem" = "CD4_Tcells",
#                                                                "CD8 Naive" = "CD8_Tcells","CD8 Mem" = "CD8_Tcells"))
seurat.obj$celltype_main <- plyr::revalue(seurat.obj$celltype, c("B Interm" = "Bcells", "B Mem" = "Bcells", "B Naive" = "Bcells",
                                                                 "CD14 Mono" = "Monocytes", "CD16 Mono" = "Monocytes",
                                                                 "Treg"= "Tcells","CD4 Naive" = "Tcells", "CD4 Mem" = "Tcells",
                                                                 "CD8 Naive" = "Tcells","CD8 Mem" = "Tcells", "MAIT"="Tcells",
                                                                 "MAIT"="Tcells", "Treg" = "Tcells", "gdT" = "Tcells", "cDC" = "Dendritic", "pDC" = "Dendritic"))
celltypes <- unique(seurat.obj$celltype_main)

DefaultAssay(seurat.obj) <- 'RNA'
seurat.obj <- NormalizeData(seurat.obj)

############################################################################
############################################################################
###                                                                      ###
###                          MULTIMODAL MARKERS                          ###
###                                                                      ###
############################################################################
############################################################################

if(!use.weights){
  DefaultAssay(seurat.obj) <- 'ATAC'
  # Extract top TFs at the single-cell level --------------------------------
  seurat.obj$size <- 1
  multimodalMarkersMethod = "weighted_t"
  # markers.summary <- FindMultimodalMarkers.SuperCell(seurat.obj = seurat.obj, group.by = 'celltype_main',
  #                                                    assay1 = "RNA", assay2 = "chromvar",
  #                                                    min.cells.feature = 0, min.cells.group = 0,
  #                                                    min.pct = 0.01, padj.cutoff = 0.05, base = 2,
  #                                                    test.use = multimodalMarkersMethod,
  #                                                    only.pos = T, fc.name1 = "avg_log2FC",
  #                                                    fc.name2 = "avg_diff")
  TFs <- ConvertMotifID(seurat.obj, id = rownames(seurat.obj[["chromvar"]]))
  motifs.subset <- rownames(seurat.obj[["chromvar"]])[TFs %in% rownames(seurat.obj[["RNA"]])]
  TFs.subset <- TFs[TFs %in% rownames(seurat.obj[["RNA"]])]
  markers.summary <- FindMultimodalMarkers.SuperCell2(seurat.obj = seurat.obj, group.by = 'celltype_main',
                                                      assay1 = "RNA", assay2 = "chromvar",
                                                      min.cells.feature = 0, min.cells.group = 0,
                                                      min.pct = 0, padj.cutoff = 1, base = 2,
                                                      logfc.threshold1 = 0, logfc.threshold2 = 0,
                                                      features.1 = TFs.subset, features.2 = motifs.subset,
                                                      test.use = multimodalMarkersMethod,
                                                      only.pos = T, fc.name1 = "avg_log2FC",
                                                      fc.name2 = "avg_diff")
  head(markers.summary[[1]])
  saveRDS(markers.summary, file = paste0(opt$outdir, "/multimodalMarkers_mainCellTypes_ttest.rds"))


  # Get TopTFs using Seurat approach ----------------------------------------

  wilcox.rna <- presto:::wilcoxauc.Seurat(
    X = seurat.obj,
    group_by = 'celltype_main',
    assay = 'data',
    seurat_assay = 'RNA'
  )
  wilcox.motifs <- presto:::wilcoxauc.Seurat(
    X = seurat.obj,
    group_by = 'celltype_main',
    assay = 'data',
    seurat_assay = 'chromvar'
  )

  motif.names <- wilcox.motifs$feature
  colnames(wilcox.rna) <- paste0("RNA.", colnames(wilcox.rna))
  colnames(wilcox.motifs) <- paste0("chromvar.", colnames(wilcox.motifs))
  wilcox.rna$gene <- wilcox.rna$RNA.feature
  DefaultAssay(seurat.obj) <- "ATAC"
  wilcox.motifs$gene <- ConvertMotifID(seurat.obj, id = motif.names)

  saveRDS(wilcox.rna, file = paste0(opt$outdir, "wilcox_rna_mainCellTypes.rds"))
  saveRDS(wilcox.motifs, file = paste0(opt$outdir, "wilcox_motifs_mainCellTypes.rds"))

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

  topTFs.wilcoxonAUC <- do.call(rbind, lapply(unique(seurat.obj$celltype_main), function(i) topTFs.original(markers_rna = wilcox.rna, markers_motifs = wilcox.motifs, celltype = i)))
  saveRDS(topTFs.wilcoxonAUC, paste0(opt$outdir,"/multimodalMarkers_mainCellTypes_wilcoxonAUC.rds"))


}else{
  DefaultAssay(seurat.obj) <- 'ATAC'

  for(multimodalMarkersMethod in c("survey_weighted_t", "weighted_t","nonWeigthed")){
    if(multimodalMarkersMethod == "nonWeigthed"){
      seurat.obj$size <- 1
      test.method = "weighted_t"
    }else{
      test.method = multimodalMarkersMethod
    }
    # markers.summary <- FindMultimodalMarkers.SuperCell(seurat.obj = seurat.obj, group.by = 'celltype_main',
    #                                                    assay1 = "RNA", assay2 = "chromvar",
    #                                                    min.cells.feature = 0, min.cells.group = 0,
    #                                                    min.pct = 0.01, padj.cutoff = 0.05, base = 2,
    #                                                    test.use = test.method,
    #                                                    only.pos = T,
    #                                                    fc.name1 = "avg_log2FC", fc.name2 = "avg_diff")
    TFs <- ConvertMotifID(seurat.obj, id = rownames(seurat.obj[["chromvar"]]))
    motifs.subset <- rownames(seurat.obj[["chromvar"]])[TFs %in% rownames(seurat.obj[["RNA"]])]
    TFs.subset <- TFs[TFs %in% rownames(seurat.obj[["RNA"]])]
    markers.summary <- FindMultimodalMarkers.SuperCell2(seurat.obj = seurat.obj, group.by = 'celltype_main',
                                                        assay1 = "RNA", assay2 = "chromvar",
                                                        min.cells.feature = 0, min.cells.group = 0,
                                                        min.pct = 0, padj.cutoff = 1, base = 2,
                                                        logfc.threshold1 = 0, logfc.threshold2 = 0,
                                                        features.1 = TFs.subset, features.2 = motifs.subset,
                                                        test.use = test.method,
                                                        only.pos = T,
                                                        fc.name1 = "avg_log2FC", fc.name2 = "avg_diff")
    head(markers.summary[[1]])
    saveRDS(markers.summary, paste0(opt$outdir,"/multimodalMarkers_mainCellTypes_", gsub("_", "",multimodalMarkersMethod), ".rds"))

  }


}
