library(SeuratData)
library(Seurat)
library(Signac)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(getopt)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds) or dataset name (from SeuratData package)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'fragmentFile', "f",1,  "character", "Fragment file path",
  "RNAcomp", "p", 1, "character", "range of components to consider for wnn analysis (eg 1:50 for RNA pca)",
  "ATACcomp", "q", 1, "character", "range of components to consider for wnn analysis (eg 2:50 for ATAC pca)",
  "k.wnn", "k", 1, "numeric", "k for the knn used in the wnn analysis",
  "RNAnormalization", "a", 1, "character", "normalisation method for RNA (logNormalize or SCTransform)", 
  "nVarGenes", "v", 1, "numeric", "number of variable genes",
  'minCutOff', "c", 1, "character", "ATAC features selection cut off (default q0)",
  "diffMap", "d", 0, "logical", "compute diffusion maps (defautl FALSE)",
  "python", "y", 1, "character", "python path with palantir installed to compute diff map"
), byrow=TRUE, ncol=5)

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code
# test
# opt <- list()
# opt$fragmentFile <- "~/Documents/multiomicsMetacells/multiome_PBMC_data/fragments_files/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
# opt$inputSeurat <- "pbmcMultiome"
# opt$outdir <- "output/correlationAnalyzis/pbmcMultiome/singlecell_analysis"
# frag.file <- opt$fragmentFile 
# opt$RNAcomp <- "1:50"
# opt$ATACcomp <- "2:50"
# opt$minCutOff <- "q0"
# opt$RNAnormalization <- "SCTransform"

if(is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "logNormalize"
}


if (is.null(opt$minCutOff)) {
  opt$minCutOff <- "q0"
}

if (is.null(opt$RNAcomp)) {
  opt$RNAcomp <- c(1:30)
} else {
  ci <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][2])
  opt$RNAcomp <- c(ci:cf)
}

if (is.null(opt$ATACcomp)) {
  opt$ATACcomp <- c(1:30)
} else {
  ci <- as.numeric(strsplit(opt$ATACcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$ATACcomp,split = ":")[[1]][2])
  opt$ATACcomp <- c(ci:cf)
}

if (is.null(opt$outdir)) {
  opt$outdir <- "./"
}

dir.create(opt$outdir,recursive = T,showWarnings = F)


if(endsWith(opt$inputSeurat,suffix = "rds")) {
  pbmc <- readRDS(opt$inputSeurat)
}else{
  data("pbmc.atac")
  
  # Now add in the ATAC-seq data
  # we'll only use peaks in standard chromosomes
  grange.counts <- StringToGRanges(rownames(pbmc.atac))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  pbmc.atac <- pbmc.atac[as.vector(grange.use), ]
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"
  
  pbmc.atac <- CreateChromatinAssay(
    counts = pbmc.atac@assays$ATAC@counts,
    genome = 'hg38',
    fragments = opt$fragmentFile,
    annotation = annotations
  )
  
  data("pbmc.rna")
  
  pbmc.rna[["ATAC"]] <- pbmc.atac
  
  remove(pbmc.atac)
  
  pbmc <- pbmc.rna
  
  remove(pbmc.rna)
}



#if present we use the cell filtering stored in the column seurat annotation
#We define coarse annotations by merging the different CD8 (resp. CD4) memory types.

if ("seurat_annotations" %in% colnames(pbmc@meta.data)) {
  pbmc <- pbmc[,pbmc$seurat_annotations != "filtered"]
  Idents(pbmc) <- "seurat_annotations"
}

addCellTypePBMC <- function(pbmc) {
  pbmc$celltype <- pbmc$seurat_annotations
  
  pbmc$celltype[grepl(pattern = "CD8 TEM",x = pbmc$celltype)] <- "CD8 Mem"
  
  pbmc$celltype[grepl(pattern = "CD4 TEM",x = pbmc$celltype)] <- "CD4 Mem"
  pbmc$celltype[grepl(pattern = "CD4 TCM",x = pbmc$celltype)] <- "CD4 Mem"
  
  pbmc$celltype[grepl(pattern = "CD8 TEM",x = pbmc$celltype)] <- "CD8 Mem"
  
  pbmc$celltype[grepl(pattern = "Intermediate B",x = pbmc$celltype)] <- "B Interm"
  pbmc$celltype[grepl(pattern = "Naive B",x = pbmc$celltype)] <- "B Naive"
  pbmc$celltype[grepl(pattern = "Memory B",x = pbmc$celltype)] <- "B Mem"
  
  Idents(pbmc) <- "celltype"
  return(pbmc)
}

pbmc <- addCellTypePBMC(pbmc)

# pbmc$coarse.annotation <- pbmc$seurat_annotations
# 
# #pbmc$coarse.annotation[grepl(pattern = "B",x = pbmc$coarse.annotation)] <- "B"
# 
# pbmc$coarse.annotation[grepl(pattern = "CD8 TEM",x = pbmc$coarse.annotation)] <- "CD8 Mem"
# 
# pbmc$coarse.annotation[grepl(pattern = "CD4 TEM",x = pbmc$coarse.annotation)] <- "CD4 Mem"
# pbmc$coarse.annotation[grepl(pattern = "CD4 TCM",x = pbmc$coarse.annotation)] <- "CD4 Mem"




# #We define a color palette for this new annotations.
# 
# color <- c("CD4 Naive"="#999999","NK"="#004949","CD8 Naive"="#009292","CD14 Mono"="#ff6db6",
#            "gdT"="#490092", "CD4 Mem"="#006ddb","cDC"="#b66dff","Treg"="#6db6ff",
#            "Intermediate B"="#b6dbff","Memory B"= "#8494FF","Naive B" = "#00A9FF",
#            "CD16 Mono"="#920000","HSPC"="#924900","CD8 Mem"="#db6d00","pDC"="#24ff24", "MAIT"="#ffff6d","Plasma"="#ffb6db")


## Analyzis of each modality separately

# Seurat scRNA-seq workflow

#As in Seurat tutorial for multimodal analyizis we use the SCTransform normalization for RNA data

if (opt$RNAnormalization == "SCTransform") {
  rnaAssay = "SCT"
  DefaultAssay(pbmc) <- "RNA"
  options(future.globals.maxSize = 8 * 1024 ^ 3) # for 50 Gb RAM
  pbmc <- SCTransform(pbmc, verbose = FALSE,conserve.memory = TRUE) %>% RunPCA() 
} else {
  rnaAssay = "RNA"
  pbmc <- NormalizeData(pbmc, verbose = FALSE) %>% FindVariableFeatures(pbmc,nFeature = opt$nVarGenes) %>% ScaleData(pbmc) %>% RunPCA()
}


# Signac scATAC-seq workflow


# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
#grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = opt$minCutOff)
pbmc <- RunSVD(pbmc)


# ## Save in h5ad for SEACells
# SeuratDisk::SaveH5Seurat(pbmc, filename =  paste0(opt$outdir,"/seurat.h5Seurat"))
# SeuratDisk::Convert(paste0(opt$outdir,"/seurat.h5Seurat"), dest =paste0(opt$outdir,"/seurat.ATAC.h5ad"),assay ="ATAC")
# SeuratDisk::Convert(paste0(opt$outdir,"/seurat.h5Seurat"), dest =paste0(opt$outdir,"/seurat.",rnaAssay,".h5ad"),assay =rnaAssay)
# if (rnaAssay != "RNA") {
# SeuratDisk::Convert(paste0(opt$outdir,"/seurat.h5Seurat"), dest =paste0(opt$outdir,"/seurat.RNA.h5ad"),assay ="RNA",overwrite = T)
# }
# system(command = paste0("rm -f ",paste0(opt$outdir,"/seurat.h5Seurat")))


## Multimodal analyzis with Seurat

pbmc <- RunUMAP(pbmc,dims = opt$RNAcomp, 
                reduction.name = 'umap.rna', 
                reduction.key = 'rnaUMAP_')

pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = opt$ATACcomp, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


# UMAP results


p1 <- DimPlot(pbmc, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("RNA")+ NoLegend()
p2 <- DimPlot(pbmc, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("ATAC")+ NoLegend()

pdf(file = paste0(opt$outdir,"/umap_single_modality.pdf"))
p1 
p2
dev.off()

pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(opt$RNAcomp, opt$ATACcomp))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",return.model = T)

#Comparison of UMAP results 

p3 <- DimPlot(pbmc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("WNN")
#p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

pdf(file = paste0(opt$outdir,"/umap_wnn_analysis.pdf"))
p3
dev.off()

pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

pdf(file = paste0(opt$outdir,"/umap_wnn_analysis_clusters.pdf"))
DimPlot(pbmc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("WNN")
dev.off()

if(opt$diffMap) {
  ## compute diffusion map for benchmarking
  library(MetacellAnalysisToolkit)
  library(reticulate)
  use_python(opt$python)
  pca_diffusion_comp <- get_diffusion_comp(sc.obj = pbmc, dims = opt$RNAcomp)
  colnames(pca_diffusion_comp) <- c(1:ncol(pca_diffusion_comp))
  pca_diffusion_comp <-  as.matrix(pca_diffusion_comp)
  pbmc[["pca_diffusion"]] <- CreateDimReducObject(embeddings = pca_diffusion_comp,assay = "RNA",key = "DM_") 
  
  lsi_diffusion_comp <- get_diffusion_comp(sc.obj = pbmc, dims = opt$ATACcomp,sc.reduction = "lsi")
  colnames(lsi_diffusion_comp) <- c(1:ncol(lsi_diffusion_comp))
  lsi_diffusion_comp <-  as.matrix(lsi_diffusion_comp)
  pbmc[["lsi_diffusion"]] <- CreateDimReducObject(embeddings = lsi_diffusion_comp,assay = "ATAC",key = "ADM_") 
}

saveRDS(pbmc, paste0(opt$outdir,"/seuratWNN.rds"))
