library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(SuperCellMultiomics)
library(anndata)
library(ggplot2)
library(getopt)

set.seed(1234)


spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'rnaPostProH5ad',  'r', 1, "character", "REQUIRED : processed (ready for moment computations) rna data with scanpy and scvelo.",
  'atacPostProH5ad',  'a', 1, "character", "REQUIRED : processed (ready for data smoothing) aggregated atac with scanpy and MultiVelo",
  'input', 'i', 1, "character", "REQUIRED : wheter to use RNA (RNA) or RNA and aggATAC (Multi) as input for metacell identification",
  'outdir',     'o', 1, "character", 'Outdir path (default ./)',
  'cellTypeLabel', "l", 1, "character", "cell type (or cluster) label in meta data",
  "gamma", "g", 1, "numeric", "gamma for metacell identificaition (default 10)",
  "scanpyScaled", "s", 0, "logical", "Specify scanpy data have been scaled",
  "singleCellUMAPs", "u", 0, "logical", "reduction method for integration, cca (default) or rpca"
), byrow=TRUE, ncol=5)


opt = getopt(spec)

#opt <- list()

# opt$gamma <- 10
# opt$rnaPostProH5ad <- 'mouse_brain_adata_rna_postpro.h5ad'
# opt$atacPostProH5ad <- 'mouse_brain_adata_atac_postpro.h5ad'
# opt$input <- "RNA"
# opt$outdir <- "./mouseSkin/fullProcessRNA"
# opt$cellTypeLabel <- "cluster"
# opt$scanpyScaled <- F
# opt$singleCellUMAPs <- T

#setwd("~/Documents/multiomicsMetacells/MultiVelo/multivelo_demo/")


# if help was asked, print a friendly message
# and exit with a non-zero error code

args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$rnaPostProH5ad)) {
  cat("integration quality control with kBET metrics on different batch labels")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
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

if (opt$scanpyScaled) {
  seuratObject <- ScaleData(seuratObject,do.scale = F,do.center = F) # already scaled in scVelo default
} else {
  seuratObject <- ScaleData(seuratObject,do.scale = T,do.center = T) # TO TEST for consistency with PCA computed with scVelo
}

seuratObject <- RunPCA(seuratObject, verbose = FALSE)

if(opt$singleCellUMAPs) {
  seuratObject <- RunUMAP(seuratObject, dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_') # optional
  
  png(paste0(opt$outdir,"/single_cells_umap.rna.png"))
  plot(DimPlot(seuratObject, reduction = "umap.rna",group.by = opt$cellTypeLabel,label = T))
  dev.off()
  
}


DefaultAssay(seuratObject) <- "aggATAC"
VariableFeatures(seuratObject) <- rownames(seuratObject)
seuratObject <- Signac::RunSVD(seuratObject, verbose = FALSE) # data are already TFIDF normalized

if (!is.null(seuratObject@meta.data$nCount_aggATAC)) {
png(paste0(opt$outdir,"/DepthCor_aggATAC_sc.png"))
Signac::DepthCor(seuratObject)
dev.off()
}

if(opt$singleCellUMAPs) {
  seuratObject <- RunUMAP(seuratObject,reduction = "lsi", dims = 1:30, reduction.name = 'umap.aggAtac', reduction.key = 'aggAtacUMAP_') # optional
  
  png(paste0(opt$outdir,"/single_cells_umap.aggAtac.png"))
  plot(DimPlot(seuratObject, reduction = "umap.aggAtac",group.by = opt$cellTypeLabel,label = T))
  dev.off()
}

if(opt$singleCellUMAPs) {
  # find weighted nearest neighbors
  seuratObject <- FindMultiModalNeighbors(seuratObject, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 1:30))
  seuratObject <- RunUMAP(seuratObject, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") # optional
  png(paste0(opt$outdir,"/single_cells_umap.wnn.png"))
  plot(DimPlot(seuratObject, reduction = "wnn.umap",group.by = opt$cellTypeLabel,label = T))
  dev.off()
}

# scanpyPCA <- object_rna$obsm$X_pca
# rownames(scanpyPCA) <- colnames(seuratObject)
# 
# seuratObject[["scanpyPCA"]] <- CreateDimReducObject(embeddings = scanpyPCA,
#                                                     key = "scPCA_",assay = "RNA")

if (opt$input == "RNA") {
  seuratMetacellObject <- SCimplify_for_Seurat(seurat = seuratObject,
                                               k.knn = 30,
                                               gamma = opt$gamma,
                                               assay = c("RNA"),  
                                               dims = list(c(1:30)),
                                               kernel = T)
}

if (opt$input == "Multi") {
  seuratMetacellObject <- SCimplify_for_Seurat(seuratObject,
                                               k.knn = 30,
                                               assay = c('RNA','aggATAC'),
                                               reduction = list("pca", "lsi"),
                                               dims = list(1:30, 1:30),
                                               graph.name = "knn",
                                               kernel = T,
                                               gamma = opt$gamma)
}

# averging UMAP coordinates
scanpyUmapCoord <- object_rna$obsm$X_umap
#scanpyUmapCoord <- cbind(scanpyUmapCoord, seuratMetacellObject@misc$membership)
colnames(scanpyUmapCoord) <- c("scUMAP_1","scUMAP_2")
centroids <- stats::aggregate(scanpyUmapCoord ~ seuratMetacellObject@misc$membership, scanpyUmapCoord, 
                              mean)
centroids <- as.matrix(centroids[,c("scUMAP_1","scUMAP_2")])
rownames(centroids) <- colnames(seuratMetacellObject)
seuratMetacellObject@meta.data[,opt$cellTypeLabel] <- factor(seuratMetacellObject@meta.data[,opt$cellTypeLabel],
                                                             levels = levels(seuratObject@meta.data[,opt$cellTypeLabel]))

seuratMetacellObject[['scUMAP']] <- CreateDimReducObject(embeddings = centroids,key = "scUMAP",assay = "RNA")

#DimPlot(seuratMetacellObject,reduction = "scUMAP",group.by = "cluster")
# Because of potential regression (RNA assay) or abscence of raw counts (aggATAC)
dataRNA <- as.matrix(supercell_GE(seuratObject[["RNA"]]@data, # contains scaled data if scanpyScaled is True
                                  groups = seuratMetacellObject@misc$membership,
                                  mode = "average"))
colnames(dataRNA) <- as.character(rownames(seuratMetacellObject@meta.data))

seuratMetacellObject[["RNA"]]@data <- dataRNA

avgAggATAC <- as.matrix(supercell_GE(seuratObject[["aggATAC"]]@data,groups = seuratMetacellObject@misc$membership,mode = "average"))
colnames(avgAggATAC) <- as.character(rownames(seuratMetacellObject@meta.data))

seuratMetacellObject[["aggATAC"]] <- CreateAssayObject(data = avgAggATAC,assay = "aggATAC")

# contracting originally computed wnn graph (ATAC connectivities of adata object)
initialWnnGraph <- igraph::graph_from_adjacency_matrix(object_atac$obsp$connectivities,mode = "undirected",weighted = T)
mcWnnGraph <- igraph::contract(initialWnnGraph, seuratMetacellObject@misc$membership)
mcWnnGraph <- igraph::simplify(mcWnnGraph, remove.loops = T, edge.attr.comb="sum")
mcAdjMatWnn <- igraph::as_adjacency_matrix(mcWnnGraph,attr = "weight")

rownames(mcAdjMatWnn) <- colnames(seuratMetacellObject)
colnames(mcAdjMatWnn) <- colnames(seuratMetacellObject)
neighbors <- FindNeighbors(as.Graph(mcAdjMatWnn),return.neighbor = T)
nn_idx <- neighbors@nn.idx
nn_dist <- neighbors@nn.dist
nn_cells <- neighbors@cell.names

# save conctracted wnn graph
write.table(nn_idx, paste0(opt$outdir,"/nn_idx.txt"), sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_dist, paste0(opt$outdir,"/nn_dist.txt"), sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_cells, paste0(opt$outdir,"/nn_cells.txt"), sep = ',', row.names = F, col.names = F, quote = F)

avgSpliced <- as.matrix(supercell_GE(seuratObject[["spliced"]]@data,groups = seuratMetacellObject@misc$membership,mode = "average"))
colnames(avgSpliced) <- as.character(rownames(seuratMetacellObject@meta.data))
seuratMetacellObject[["spliced"]] <- CreateAssayObject(data = avgSpliced,assay = "spliced")

avgUnspliced <- as.matrix(supercell_GE(seuratObject[["unspliced"]]@data,groups = seuratMetacellObject@misc$membership,mode = "average"))
colnames(avgUnspliced) <- as.character(rownames(seuratMetacellObject@meta.data))

#seuratMetacellObject[["unspliced"]] <- CreateAssayObject(data = avgUnspliced,assay = "unspliced")

avgAggATAC <- as.matrix(supercell_GE(seuratObject[["aggATAC"]]@data,groups = seuratMetacellObject@misc$membership,mode = "average"))
colnames(avgAggATAC) <- as.character(rownames(seuratMetacellObject@meta.data))

#seuratMetacellObject[["aggATAC"]] <- CreateAssayObject(data = avgAggATAC,assay = "aggATAC")

filePrefix <- paste0(opt$outdir,"/adata_mc")

write_h5ad(AnnData(X = t(seuratMetacellObject[["aggATAC"]]@data),
                   obs = seuratMetacellObject@meta.data,
                   obsp = list("connectivities" = mcAdjMatWnn)),
           filename = paste0(filePrefix,"_atac.h5ad"))

write_h5ad(AnnData(X = t(seuratMetacellObject[["RNA"]]@data),
                   obs = seuratMetacellObject@meta.data,
                   obsm = list('X_umap'= seuratMetacellObject@reductions$scUMAP@cell.embeddings),
                   layers = list("spliced" = t(avgSpliced),"unspliced" = t(avgUnspliced))),
           filename = paste0(filePrefix,"_rna.h5ad"))


