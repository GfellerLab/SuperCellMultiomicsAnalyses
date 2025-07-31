library(Seurat)
library(dplyr)
library(getopt)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "sampleName",  "s", 1, "character", "sample name",
  "RNAcomp", "p", 1, "character", "range of components to consider for wnn analysis (eg 1:50 for RNA pca)",
  "ADTcomp", "q", 1, "character", "range of components to consider for wnn analysis (eg 2:50 for ATAC pca)",
  "diffMap", "d", 0, "logical", "compute diffusion maps (defautl FALSE)",
  "python", "y", 1, "character", "python path with palantir installed to compute diff map",
  "RNAnormalization", "r", 1, "character", "normalisation method for RNA (logNormalize or SCTransform)"

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

if(is.null(opt$diffMap)) {
  opt$diffMap <- F

}

if(is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "logNormalize"

}


if (is.null(opt$RNAcomp)) {
  opt$RNAcomp <- c(1:30)
} else {
  ci <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][2])
  opt$RNAcomp <- c(ci:cf)
}

if (is.null(opt$ADTcomp)) {
  opt$ADTcomp <- c(1:30)
} else {
  ci <- as.numeric(strsplit(opt$ADTcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$ADTcomp,split = ":")[[1]][2])
  opt$ADTcomp <- c(ci:cf)
}

if (is.null(opt$outdir)) {
  opt$outdir <- "./"
}

print(opt)

dir.create(opt$outdir,recursive = T,showWarnings = F)
options(future.globals.maxSize = 2000 * 1024^2) 

pbmc <- readRDS(opt$inputSeurat)
pbmc$orig.ident <- paste0(pbmc$donor,"_",pbmc$time)
pbmc <- pbmc[,pbmc$orig.ident == opt$sampleName]
pbmc <- pbmc[,pbmc$celltype.l2 != "Doublet"]
#pbmc <- pbmc[,pbmc$time == 0]  # first test on reference samples
pbmc$celltype.l1.5 <- pbmc$celltype.l1
pbmc$celltype.l1.5[pbmc$celltype.l1 == "other"] <- pbmc$celltype.l2[pbmc$celltype.l1 == "other"]

pbmc$celltype.l1.5[pbmc$celltype.l1 == "other T"] <- pbmc$celltype.l2[pbmc$celltype.l1 == "other T"]

gc()

DefaultAssay(pbmc) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the
VariableFeatures(pbmc) <- rownames(pbmc[["ADT"]])
pbmc <- NormalizeData(pbmc, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca')

# DefaultAssay(pbmc) <- 'RNA'
#

DefaultAssay(pbmc) <- "RNA"
if (opt$RNAnormalization == "SCT") {
  rnaAssay = "SCT"
  pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA()
} else {
  rnaAssay = "RNA"
  pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",assay = "RNA") %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
}

## Save in h5ad for SEACells
# SeuratDisk::SaveH5Seurat(pbmc, filename =  paste0(opt$outdir,"/seurat.h5Seurat"),overwrite = T)
# SeuratDisk::Convert(paste0(opt$outdir,"/seurat.h5Seurat"), dest =paste0(opt$outdir,"/seurat.ADT.h5ad"),assay ="ADT",overwrite = T)
# SeuratDisk::Convert(paste0(opt$outdir,"/seurat.h5Seurat"), dest =paste0(opt$outdir,"/seurat.",rnaAssay,".h5ad"),assay =rnaAssay,overwrite = T)
# if (rnaAssay != "RNA") {
# SeuratDisk::Convert(paste0(opt$outdir,"/seurat.h5Seurat"), dest =paste0(opt$outdir,"/seurat.RNA.h5ad"),assay ="RNA",overwrite = T)
# }
# system(command = paste0("rm -f ",paste0(opt$outdir,"/seurat.h5Seurat")))

pbmc <- FindMultiModalNeighbors(
  pbmc,
  reduction.list = list("pca", "apca"),
  dims.list = list(opt$RNAcomp, opt$ADTcomp),
  modality.weight.name = "RNA.weight",
  return.intermediate = T
)
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

pdf(file = paste0(opt$outdir,"/pbmc_Cite_wnn_umap.pdf"))
DimPlot(pbmc,reduction = "wnn.umap", group.by = 'celltype.l1')
DimPlot(pbmc,reduction = "wnn.umap", group.by = 'celltype.l2')
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

  apca_diffusion_comp <- get_diffusion_comp(sc.obj = pbmc, dims = opt$ADTcomp,sc.reduction = "apca")
  colnames(apca_diffusion_comp) <- c(1:ncol(apca_diffusion_comp))
  apca_diffusion_comp <-  as.matrix(apca_diffusion_comp)
  pbmc[["apca_diffusion"]] <- CreateDimReducObject(embeddings = apca_diffusion_comp,assay = "ADT",key = "ADM_")
}

saveRDS(pbmc, paste0(opt$outdir,"/seuratWNN.rds"))
