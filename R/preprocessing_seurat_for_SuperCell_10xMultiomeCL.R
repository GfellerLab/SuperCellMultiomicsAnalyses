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



if (opt$RNAnormalization == "SCTransform") {
  rnaAssay = "SCT"
  DefaultAssay(pbmc) <- "RNA"
  options(future.globals.maxSize = 8 * 1024 ^ 3) # for 50 Gb RAM

  pbmc <- SCTransform(pbmc, verbose = FALSE,conserve.memory = T) %>% RunPCA()
} else {
  rnaAssay = "RNA"
  pbmc <- NormalizeData(pbmc, verbose = FALSE) %>% FindVariableFeatures(pbmc,nFeature = opt$nVarGenes) %>% ScaleData(pbmc) %>% RunPCA()
}


DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = opt$minCutOff)
pbmc <- RunSVD(pbmc)




saveRDS(pbmc, paste0(opt$outdir,"/seurat_multimodal.rds"))
