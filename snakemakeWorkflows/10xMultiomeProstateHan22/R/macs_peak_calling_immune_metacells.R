library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(EnsDb.Mmusculus.v79)
library(Signac)
library(S4Vectors)
library(patchwork)
library(getopt)


options(future.globals.maxSize= 8000*1024^2)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'input',     'i',1, "character", 'input RDS',
  'outdir',     'o',1, "character", 'Outdir path (default ./)'
), byrow=TRUE, ncol=5)

opt = getopt(spec)

print(opt)

combined.metacells <- readRDS(opt$input)



DepthCor(combined.metacells)



DefaultAssay(combined.metacells) <- "ATAC"
peaks <- CallPeaks(
  object = combined.metacells,
  #macs2.path = "/home/leonard/bin/miniconda3/bin/macs3",
  group.by = "major_type"
)

peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

saveRDS(peaks,paste0(opt$outdir,"/peaks.rds"))