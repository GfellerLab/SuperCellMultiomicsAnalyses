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
  'input',  'i', 1, "character", "REQUIRED : seurat object (.rds) or adata object (h5ad)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'fragmentFile', "f",1,  "character", "Fragment file path (Z rep sep by+)",
  "k.wnn", "k", 1, "numeric", "k for the knn used in the wnn analysis",
  "nVarGenes", "v", 1, "numeric", "number of variable genes",
  'minCutOff', "c", 1, "character", "ATAC features selection cut off (default q0)"

), byrow=TRUE, ncol=5)

opt = getopt(spec)


if(is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "logNormalize"
}

if(!is.null(opt$fragmentFile)) {
  opt$fragmentFile <- strsplit(opt$fragmentFile,"\\+")[[1]]
}

if (is.null(opt$minCutOff)) {
  opt$minCutOff <- "q0"
}


if (is.null(opt$outdir)) {
  opt$outdir <- "./"
}

dir.create(opt$outdir,recursive = T,showWarnings = F)

print(opt)

if(endsWith(opt$input,suffix = "rds")) {
  hspc <- readRDS(opt$input)
}

if(endsWith(opt$input,suffix = ".h5ad")){
  atac <- anndata::read_h5ad(opt$input)

  atac.count <- Matrix::t(atac$X)

  rep1.ori <- colnames(atac.count)[grepl(x = colnames(atac.count),pattern = "rep1")]
  rep1 <- sub(x= rep1.ori,"cd34_multiome_rep1#",replacement = "")

  rep2.ori <- colnames(atac.count)[grepl(x = colnames(atac.count),pattern = "rep2")]
  rep2 <- sub(x= rep2.ori,"cd34_multiome_rep2#",replacement = "")
  names(rep2) <- rep2.ori
  names(rep1) <- rep1.ori

  grange.counts <- StringToGRanges(rownames(atac.count), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac.count <- atac.count[as.vector(grange.use), ]
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"

  chrom_assay <- CreateChromatinAssay(
    counts = atac.count,
    sep = c(":", "-"),
    genome = 'hg38',
    annotation = annotations
  )

  hspc <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    meta.data = atac$obs
  )
  message("atac object created...")

  remove(atac.count)
  gc()

}




# Signac scATAC-seq workflow


# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
#grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
DefaultAssay(hspc) <- "ATAC"
hspc <- RunTFIDF(hspc)
hspc <- FindTopFeatures(hspc, min.cutoff = opt$minCutOff)
hspc <- RunSVD(hspc)



## Save in h5ad for SEACells
adata <- anndata::AnnData(X = Matrix::t(GetAssayData(object = hspc,slot = "counts",assay = "ATAC")),
                          obs = hspc@meta.data,
                          #raw = adata.raw,
                          obsm = list("X_lsi" = hspc[["lsi"]]@cell.embeddings))

anndata::write_h5ad(adata,paste0(opt$outdir,"/seurat.ATAC.h5ad"))
