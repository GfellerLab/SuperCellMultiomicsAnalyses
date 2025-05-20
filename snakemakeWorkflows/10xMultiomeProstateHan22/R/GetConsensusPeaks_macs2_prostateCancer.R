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
  'indir',        'i', 1, "character", "inputDir the peak folder (one per smp)",
  'outdir',     'o', 1, "character", 'Outdir path (default ./)'
), byrow=TRUE, ncol=5)

opt = getopt(spec)


peakFiles <- list.files(path =opt$indir,pattern = ".rds",recursive = T,full.names = T)

peak.list <- lapply(X =peakFiles,FUN = function(x) {
  x <- readRDS(x)
})

combined.peaks <- GenomicRanges::reduce(x=c(peak.list[[1]],peak.list[[2]]))
for (p in peak.list[2:length(peak.list)]) {
  combined.peaks <- GenomicRanges::reduce(x=c(combined.peaks,p))
}

peakwidths <- GenomicRanges::width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

write.csv(data.frame(combined.peaks),paste0(opt$outdir,"/peakset_macs2.csv"))
