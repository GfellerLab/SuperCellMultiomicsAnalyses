library(SeuratData)
library(Seurat)
library(Signac)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SuperCellMultiomics)
library(getopt)
library(doParallel)
library(chromVAR) #
library(JASPAR2020) #
library(TFBSTools) #
library(motifmatchr) #
library(ggplot2)


spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'pwm', "p", 1, "character", "pwm motif matrix for chromvar",
), byrow=TRUE, ncol=5)


opt = getopt(spec)

opt <- list()
#setwd("~/work/SuperCellMultiomicsAnalyses/")
opt$singleCellSeurat <- "output/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
opt$pwm <- "input/JASPAR_2024_human_motifs/pwm.rds"
opt$outdir <- 'output/pbmcMultiome/singlecells_analysis/'




print(opt)

seurat.sc <- readRDS(opt$singleCellSeurat)





DefaultAssay(seurat.sc) <- "ATAC"
peaks <- CallPeaks(
  object = seurat.sc,
  #macs2.path = "/home/leonard/bin/miniconda3/envs/MACS_env/bin/macs2",
  group.by = "wsnn_res.0.8" # need to be computed before
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


seurat.sc[['peaks']] <- chrom_assay
saveRDS(paste0(opt$outdir,"/seuratWNN_with_macs2_peaks.rds"))


