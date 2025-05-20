library(Seurat)
library(dplyr)
library(getopt)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds or name for laoding with SeuratData)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "RNAcomp", "p", 1, "character", "range of components to consider for wnn analysis (eg 1:30 for RNA pca)",
  "ADTcomp", "q", 1, "character", "range of components to consider for wnn analysis (eg 1:18 for ADT pca)",
  "RNAnormalization", "r", 1, "character", "normalisation method for RNA (logNormalize or SCTransform)",
  "diffMap", "d", 0, "logical", "compute diffusion maps (defautl FALSE)",
  "python", "y", 1, "character", "python path with palantir installed to compute diff map"

), byrow=TRUE, ncol=5)

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code
# test
# opt <- list()
# opt$fragmentFile <- "~/Documents/multiomicsMetacells/multiome_bm_data/fragments_files/bm_granulocyte_sorted_10k_atac_fragments.tsv.gz"
# opt$inputSeurat <- "bmMultiome"
# opt$outdir <- "output/correlationAnalyzis/bmMultiome/singlecell_analysis"
# frag.file <- opt$fragmentFile 
# opt$RNAcomp <- "1:50"
# opt$ATACcomp <- "2:50"
# opt$minCutOff <- "q0"
# opt$RNAnormalization <- "SCTransform"

if(is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "logNormalize"
  
}

if(is.null(opt$diffMap)) {
  opt$diffMap <- F
  
}



if (is.null(opt$RNAcomp)) {
  opt$RNAcomp <- c(1:30)
} else {
  ci <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][2])
  opt$RNAcomp <- c(ci:cf)
}

if (is.null(opt$ADTcomp)) {
  opt$ADTcomp <- c(1:18)
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

if (endsWith(opt$inputSeurat,'.rds')) {
bm <- readRDS(opt$inputSeurat)
} else {
  bm <- SeuratData::LoadData(ds = opt$inputSeurat)

}

DefaultAssay(bm) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

# DefaultAssay(bm) <- 'RNA'
# 

DefaultAssay(bm) <- "RNA"
if (opt$RNAnormalization == "SCT") {
  rnaAssay = "SCT"
  bm <- SCTransform(bm, verbose = FALSE) %>% RunPCA() 
} else {
  rnaAssay = "RNA"
  bm <- NormalizeData(bm,normalization.method = "LogNormalize",assay = "RNA") %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
}


bm <- FindMultiModalNeighbors(
  bm, 
  reduction.list = list("pca", "apca"), 
  dims.list = list(opt$RNAcomp, opt$ADTcomp), 
  modality.weight.name = "RNA.weight",
  return.intermediate = T
)
bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

pdf(file = paste0(opt$outdir,"/bm_Cite_wnn_umap.pdf")) 
DimPlot(bm,reduction = "wnn.umap", group.by = 'celltype.l1')
DimPlot(bm,reduction = "wnn.umap", group.by = 'celltype.l2')
dev.off()

if(opt$diffMap) {
## compute diffusion map for benchmarking
  library(MetacellAnalysisToolkit)
  library(reticulate)
  use_python(opt$python)
  pca_diffusion_comp <- get_diffusion_comp(sc.obj = bm, dims = opt$RNAcomp)
  colnames(pca_diffusion_comp) <- c(1:ncol(pca_diffusion_comp))
  pca_diffusion_comp <-  as.matrix(pca_diffusion_comp)
  bm[["pca_diffusion"]] <- CreateDimReducObject(embeddings = pca_diffusion_comp,assay = "RNA",key = "DM_") 
  
  apca_diffusion_comp <- get_diffusion_comp(sc.obj = bm, dims = opt$ADTcomp,sc.reduction = "apca")
  colnames(apca_diffusion_comp) <- c(1:ncol(apca_diffusion_comp))
  apca_diffusion_comp <-  as.matrix(apca_diffusion_comp)
  bm[["apca_diffusion"]] <- CreateDimReducObject(embeddings = apca_diffusion_comp,assay = "ADT",key = "ADM_") 
}

saveRDS(bm, paste0(opt$outdir,"/seuratWNN.rds"))

