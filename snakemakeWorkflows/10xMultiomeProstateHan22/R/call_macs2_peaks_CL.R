{library(getopt)
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(EnsDb.Mmusculus.v79)
  library(Signac)
  library(S4Vectors)
  library(patchwork)
  set.seed(1234)}


spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputH5',  'i', 1, "character", "REQUIRED : seurat object (.rds) or dataset name (from SeuratData package)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'fragmentFile', "f",1,  "character", "Fragment file path"
), byrow=TRUE, ncol=5)

opt = getopt(spec)


#if (is.null(opt$outdir)) {opt$outdir = "output/10xMultiomeProstateHan22/"}

print(opt)

smp <- strsplit(opt$inputH5,"/")[[1]][length(strsplit(opt$inputH5,split = "/")[[1]])-1]
inputdata.10x <- Read10X_h5(opt$inputH5)
grange.counts<- StringToGRanges(rownames(inputdata.10x$Peaks), sep = c(":", "-"))
grange.use<- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts<- inputdata.10x$Peaks[as.vector(grange.use), ]
annotations<- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
fragment <- opt$fragmentFile
chrom_assay<- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = fragment,
  min.cells = 10,
  annotation = annotations)

peaks <- CallPeaks(chrom_assay)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

saveRDS(peaks,paste0(opt$outdir,"/macs2_peaks.rds"))












