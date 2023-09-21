library(chromVAR) #
library(JASPAR2020) #
library(TFBSTools) #
library(motifmatchr) #
library(BSgenome.Hsapiens.UCSC.hg38) #
library(Signac) #
library(EnsDb.Hsapiens.v86) #
library(SeuratData)
library(Seurat)
library(Signac)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SuperCellMultiomics)
library(getopt)
library(doParallel)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds) or dataset name (from SeuratData package)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'nWorkers', 'w', 1, 'numeric'," 'number of cores to use"
), byrow=TRUE, ncol=5)

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code
# test
opt <- list()
opt$inputSeurat <- "output/correlationAnalyzis/pbmcMultiome/metacells/g200/seurat.multiome.mc.activity.rds"
opt$outdir <- "output/correlationAnalyzis/pbmcMultiome/metacells/g200/"

if(is.null(opt$nWorkers)){
  opt$nWorkers = parallel::detectCores()-4
} 

if(is.null(opt$outdir)) {
  opt$outdir = "./"
}



print(opt)

dir.create(opt$outdir,recursive = T,showWarnings = F)


seurat <- readRDS(opt$inputSeurat)

## ChromVar analysis


library(BiocParallel)
register(MulticoreParam(opt$nWorkers))

DefaultAssay(seurat) <- "ATAC"


# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(seurat), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
seurat <- SetAssayData(seurat, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes 
seurat <- RunChromVAR(
  object = seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(seurat,file = paste0(opt$outdir,"/seurat.mc.chromVar_activities.rds"))
