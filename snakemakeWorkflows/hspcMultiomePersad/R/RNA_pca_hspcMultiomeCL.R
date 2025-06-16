library(reticulate)
use_python("/opt/conda/envs/MetacellAnalysisToolkit/bin/python", required = TRUE)
library(SeuratData)
library(Seurat)
library(Signac)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(getopt)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'input',  'i', 1, "character", "REQUIRED : seurat object (.rds) or anndata object (.h5ad)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'fragmentFile', "f",1,  "character", "Fragment file path",
  "k.wnn", "k", 1, "numeric", "k for the knn used in the wnn analysis",
  "RNAnormalization", "a", 1, "character", "normalisation method for RNA (logNormalize or SCTransform)",
  "nVarGenes", "v", 1, "numeric", "number of variable genes",
  'minCutOff', "c", 1, "character", "ATAC features selection cut off (default q0)"

), byrow=TRUE, ncol=5)

opt = getopt(spec)


if(is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "logNormalize"

}

if (is.null(opt$minCutOff)) {
  opt$minCutOff <- "q0"
}



if (is.null(opt$outdir)) {
  opt$outdir <- "./"
}

dir.create(opt$outdir,recursive = T,showWarnings = F)


if(endsWith(opt$input,suffix = "h5ad")) {
  rna <- anndata::read_h5ad(opt$input)

  rna.count <- Matrix::t(rna$raw$X)
  rownames(rna.count) <- rna$var_names
  colnames(rna.count) <- rna$obs_names

  hspc <- CreateSeuratObject(counts = rna.count,meta.data = rna$obs)
}

if(endsWith(opt$input,suffix = "rds")) {
  hspc <- readRDS(opt$input)
}



if (opt$RNAnormalization == "SCTransform") {
  rnaAssay = "SCT"
  DefaultAssay(hspc) <- "RNA"
  hspc <- SCTransform(hspc, verbose = FALSE,conserve.memory = T) %>% RunPCA()
} else {
  rnaAssay = "RNA"
  hspc <- NormalizeData(hspc, verbose = FALSE) %>% FindVariableFeatures(nFeature = opt$nVarGenes) %>% ScaleData() %>% RunPCA()
}


adata <- anndata::AnnData(X = Matrix::t(GetAssayData(object = hspc,slot = "counts",assay = "RNA")),
                 obs = hspc@meta.data,
                 #raw = adata.raw,
                 obsm = list("X_pca" = hspc[["pca"]]@cell.embeddings))

anndata::write_h5ad(adata,paste0(opt$outdir,"/seurat.RNA.h5ad"))
