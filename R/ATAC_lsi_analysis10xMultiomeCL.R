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
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds) or dataset name (from SeuratData package)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'fragmentFile', "f",1,  "character", "Fragment file path",
  "k.wnn", "k", 1, "numeric", "k for the knn used in the wnn analysis",
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

print(opt)

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

  chrom_assay <- CreateChromatinAssay(
    counts = pbmc.atac@assays$ATAC@counts,
    genome = 'hg38',
    fragments = opt$fragmentFile,
    annotation = annotations
  )
  pbmc <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    meta.data = pbmc.atac@meta.data
  )
  message("atac object created...")

  remove(pbmc.atac)
  gc()

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

#if present we use the cell filtering stored in the column seurat annotation
#We define coarse annotations by merging the different CD8 (resp. CD4) memory types.

if ("seurat_annotations" %in% colnames(pbmc@meta.data)) {
  pbmc <- pbmc[,pbmc$seurat_annotations != "filtered"]
  Idents(pbmc) <- "seurat_annotations"
}



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




# Signac scATAC-seq workflow


# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
#grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = opt$minCutOff)
pbmc <- RunSVD(pbmc)

adata <- anndata::AnnData(X = Matrix::t(GetAssayData(object = pbmc,slot = "counts",assay = "ATAC")),
                          obs = pbmc@meta.data,
                          #raw = adata.raw,
                          obsm = list("X_lsi" = pbmc[["lsi"]]@cell.embeddings))

anndata::write_h5ad(adata,paste0(opt$outdir,"/seurat.ATAC.h5ad"))

## Save in h5ad for SEACells
# SeuratDisk::SaveH5Seurat(pbmc, filename =  paste0(opt$outdir,"/seurat.ATAC.h5Seurat"))
# SeuratDisk::Convert(paste0(opt$outdir,"/seurat.ATAC.h5Seurat"), dest =paste0(opt$outdir,"/seurat.ATAC.h5ad"),assay ="ATAC")
# system(command = paste0("rm -f ",paste0(opt$outdir,"/seurat.ATAC.h5Seurat")))
