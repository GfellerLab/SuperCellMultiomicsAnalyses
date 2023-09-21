library(Seurat)
library(Signac)
#library(SeuratObject)
#library(SeuratDisk)
library(SuperCellMultiomics)
library(anndata)
library(ggplot2)
library(getopt)

set.seed(2023)


spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'rnaPostProH5ad',  'r', 1, "character", "REQUIRED : processed (ready for moment computations) rna data with scanpy and scvelo.",
  'atacPostProH5ad',  'a', 1, "character", "REQUIRED : processed (ready for data smoothing) aggregated atac with scanpy and MultiVelo",
  'rawRNA',  't', 1, "character", "REQUIRED : raw rna data to compute wnn for atac smoothing.",
  'rawATAC',  'b', 1, "character", "REQUIRED : raw atac data wto compute wnn for atac smoothing",
  'outdir',     'o', 1, "character", 'Outdir path (default ./)',
  'cellTypeLabel', "l", 1, "character", "cell type (or cluster) label in meta data",
  "gamma", "g", 1, "numeric", "gamma for metacell identificaition (default 10)",
  'knn', "k", 1, "numeric", "k for wnn graph on subsampled data (default 50)",
  "scanpyScaled", "s", 0, "logical", "Specify scanpy data have been scaled for cell cycle correction",
  "singleCellUMAPs", "u", 0, "logical", "reduction method for integration, cca (default) or rpca"
), byrow=TRUE, ncol=5)


opt = getopt(spec)

# opt <- list()
# # 
# setwd("~/work/SuperCellMultiomicsAnalyses/")
# opt$gamma <- 10
# opt$rnaPostProH5ad <- 'input/multiomicVelocity/humanHSPC/adata_postpro_rna.h5ad'
# opt$rawRNA <- "input/multiomicVelocity/humanHSPC/adata_raw_rna.h5ad"
# opt$rawATAC <- "input/multiomicVelocity/humanHSPC/adata_raw_atac.h5ad"
# opt$atacPostProH5ad <- 'input/multiomicVelocity/humanHSPC/adata_postpro_atac.h5ad'
# opt$outdir <- "input/multiomicVelocity/humanHSPC/Subsampling/"
# opt$cellTypeLabel <- "leiden"
# opt$scanpyScaled <- T
# opt$singleCellUMAPs <- T
# opt$rawRNA
# 
# setwd("~/Documents/multiomicsMetacells/MultiVelo/multivelo_demo/")


# if help was asked, print a friendly message
# and exit with a non-zero error code

args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$rnaPostProH5ad)) {
  cat("integration quality control with kBET metrics on different batch labels")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if (is.null(opt$knn)) {
  opt$knn <- 50
}

if (is.null(opt$outdir)) {
  opt$outdir = "./"
}

if (is.null(opt$scanpyScaled)) {
  opt$scanpyScaled = F
}

if (is.null(opt$singleCellUMAPs)) {
  opt$singleCellUMAPs = F
}

print(opt)

dir.create(opt$outdir,recursive = T)
# Load 10X raw data

object_rna <- read_h5ad(opt$rnaPostProH5ad)
object_atac <- read_h5ad(opt$atacPostProH5ad)

print("scanpy objects loaded... ")

if (!is.null(object_rna$layers["matrix"])) {
  seuratObject <- CreateSeuratObject(
    counts = Matrix::t(object_rna$layers["matrix"]),
    meta = object_rna$obs,
    assay = "RNA") 
} else {
  seuratObject <- CreateSeuratObject(
    counts = Matrix::t(object_rna$X),
    meta = object_rna$obs,
    assay = "RNA") 
}

seuratObject[["RNA"]]@data <- as.matrix(Matrix::t(object_rna$X))
seuratObject[["spliced"]] <- CreateAssayObject(data = as.matrix(Matrix::t(object_rna$layers["spliced"])),
                                               assay = "spliced")
seuratObject[["unspliced"]] <- CreateAssayObject(data = as.matrix(Matrix::t(object_rna$layers["unspliced"])),
                                                 assay = "unspliced")


seuratObject[["aggATAC"]] <- CreateAssayObject(
  data = as.matrix(Matrix::t(object_atac$X)),
  assay = "aggATAC")

seuratObject@meta.data$nCount_aggATAC <- object_atac$obs$n_counts 

DefaultAssay(seuratObject) <- "RNA"
VariableFeatures(seuratObject) <- rownames(seuratObject)

# adding original UMAP coordinates
scanpyUmapCoord <- object_rna$obsm$X_umap
#scanpyUmapCoord <- cbind(scanpyUmapCoord, seuratMetacellObject@misc$membership)
colnames(scanpyUmapCoord) <- c("scUMAP_1","scUMAP_2")
rownames(scanpyUmapCoord) <- colnames(seuratObject)
seuratObject[['scUMAP']] <- CreateDimReducObject(embeddings = scanpyUmapCoord,key = "scUMAP",assay = "RNA")


sampledCells <- sample(colnames(seuratObject),size = length(colnames(seuratObject))/opt$gamma)

seuratObject <- subset(seuratObject,cells = sampledCells)



if(opt$singleCellUMAPs) {
  png(paste0(opt$outdir,"/sub_single_cells_umap.png"))
  plot(DimPlot(seuratObject, reduction = "scUMAP",group.by = opt$cellTypeLabel,label = T))
  dev.off()
}

filePrefix <- paste0(opt$outdir,"/adata_mc")

write_h5ad(AnnData(X = t(seuratObject[["aggATAC"]]@data),
                   obs = seuratObject@meta.data),
           #obsp = list("connectivities" = seuratObject@graphs$wknn)),
           filename = paste0(filePrefix,"_atac.h5ad"))

write_h5ad(AnnData(X = t(seuratObject[["RNA"]]@data),
                   obs = seuratObject@meta.data,
                   obsm = list('X_umap'= seuratObject@reductions$scUMAP@cell.embeddings),
                   layers = list("spliced" = t(seuratObject[["spliced"]]@data),"unspliced" = t(seuratObject[["unspliced"]]@data))),
           filename = paste0(filePrefix,"_rna.h5ad"))

# Recomputing wnn graph for smoothing aggregated atac data


object_rna <- read_h5ad(opt$rawRNA)[sampledCells,]
object_atac <- read_h5ad(opt$rawATAC)[sampledCells,]

print("scanpy objects loaded... ")

seuratObject <- CreateSeuratObject(
  counts = Matrix::t(object_rna$X),
  meta = object_rna$obs,
  assay = "RNA")
seuratObject <- NormalizeData(seuratObject)
seuratObject <- FindVariableFeatures(seuratObject)

if (opt$scanpyScaled) {
  seuratObject <- CellCycleScoring(seuratObject,g2m.features = cc.genes$g2m.genes,
                                   s.features = cc.genes$s.genes)
  seuratObject <- ScaleData(seuratObject,vars.to.regress = c("G2M.Score","S.Score")) 
  
} else {
  seuratObject <- ScaleData(seuratObject,do.scale = F) 
}
# scaled only when cell cycle is corrected for consistency with scVelo (optionally, use SCTransform)
seuratObject <- RunPCA(seuratObject, verbose = FALSE)

if (opt$singleCellUMAPs) {
seuratObject <- RunUMAP(seuratObject, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_') # optional
DimPlot(seuratObject,group.by = opt$cellTypeLabel,label = T)
}


# preprocess ATAC
seuratObject[["ATAC"]] <- CreateAssayObject(counts = Matrix::t(object_rna$X), min.cells = 1)
DefaultAssay(seuratObject) <- "ATAC"
seuratObject <- RunTFIDF(seuratObject)
seuratObject <- FindTopFeatures(seuratObject, min.cutoff = 'q0')
seuratObject <- RunSVD(seuratObject)

if (opt$singleCellUMAPs) {
  seuratObject <- RunUMAP(seuratObject, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_") # optional
  DimPlot(seuratObject,group.by = opt$cellTypeLabel,label = T)
}

# find weighted nearest neighbors
seuratObject <- FindMultiModalNeighbors(seuratObject, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50),k.nn = opt$knn)
seuratObject <- RunUMAP(seuratObject, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") # optional
if (opt$singleCellUMAPs) {
  png(paste0(opt$outdir,"/sub_single_cells_wnn_umap.png"))
  seuratObject <- RunUMAP(seuratObject, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") # optional
  plot(DimPlot(seuratObject,group.by = opt$cellTypeLabel,label = T,reduction = "wnn.umap"))
  dev.off()
}


# save new wnn graph
write.table(seuratObject@neighbors$weighted.nn@nn.idx, paste0(opt$outdir,"/nn_idx.txt"), sep = ',', row.names = F, col.names = F, quote = F)
write.table(seuratObject@neighbors$weighted.nn@nn.dist, paste0(opt$outdir,"/nn_dist.txt"), sep = ',', row.names = F, col.names = F, quote = F)
write.table(seuratObject@neighbors$weighted.nn@cell.names, paste0(opt$outdir,"/nn_cells.txt"), sep = ',', row.names = F, col.names = F, quote = F)
