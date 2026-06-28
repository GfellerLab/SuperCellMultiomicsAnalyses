library(SeuratData)
library(Seurat)
library(Signac)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(getopt)
library(doParallel)
library(ggplot2)


spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'singleCellSeurat',     'i', 1, "character", "input path for seurat sc object",
  'outdir',     'o',1, "character", 'Outdir path (default ./)'
), byrow=TRUE, ncol=5)


opt = getopt(spec)

print(opt)

seurat.sc <- readRDS(opt$singleCellSeurat)


DefaultAssay(seurat.sc) <- "ATAC"
peaks <- CallPeaks(
  object = seurat.sc,group.by = "seurat_annotations",
  additional.args = "--max-gap 50"
  #macs2.path = "/home/leonard/bin/miniconda3/envs/MACS_env/bin/macs2",
  #group.by = "wsnn_res.0.8" # need to be computed before
)

peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

peak.counts <- FeatureMatrix(
  fragments = Fragments(seurat.sc),
  features = peaks,
  cells = colnames(seurat.sc)
)


annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

chrom_assay<- CreateChromatinAssay(
  counts = peak.counts,
  #sep = c(":", "-"),
  genome = 'hg38',
  min.cells = 1, # peaks have been identified after a first analysis
  annotation=annotations,
  fragments = Fragments(seurat.sc)
)

DefaultAssay(seurat.sc) <- "RNA"
seurat.sc[["ATAC"]] <- NULL
seurat.sc[['ATAC']] <- chrom_assay
saveRDS(seurat.sc,paste0(opt$outdir,"/pbmcMultiome.rds"))
