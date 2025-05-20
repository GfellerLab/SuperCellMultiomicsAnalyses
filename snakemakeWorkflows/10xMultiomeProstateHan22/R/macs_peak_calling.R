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
library(BSgenome.Mmusculus.UCSC.mm10)


options(future.globals.maxSize= 8000*1024^2)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'input',     'i',1, "character", 'input RDS',
  'meta',       'm', 1, "character", "metadata table containing major_type annotation",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'output',      'u', 1, 'character',"type of output (seurat object default or peaks)",
  "motifs",    "p", 1, "character", "path to pfm matrix (.rds) or JASPAR version, need to be filled in order to run a chromvar analysis",
  "cores",      'c', 1, "character", "number of cores",
  "assay",      "a", 1, "character", "which chromatine assay to use for chromvar analysis (default peaks for macs2 peaks)"
  
), byrow=TRUE, ncol=5)

opt = getopt(spec)

if (is.null(opt$output)) {
  opt$output <- "seurat"
}

print(opt)

combined.metacells <- readRDS(opt$input)
allMeta <- read.csv(opt$meta)
combined.metacells$major_type <- allMeta[colnames(combined.metacells),"major_type"]


DepthCor(combined.metacells)



DefaultAssay(combined.metacells) <- "ATAC"
peaks <- CallPeaks(
  object = combined.metacells,
  #macs2.path = "/home/leonard/bin/miniconda3/bin/macs3",
  group.by = "major_type"
)

peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

if (opt$output == "seurat") {
  
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(combined.metacells),
    features = peaks,
    cells = colnames(combined.metacells)
  )
  
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "mm10"
  
  combined.metacells[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    genome = 'mm10',
    min.cells = 1,
    fragments = Fragments(combined.metacells),
    annotation = annotations
  )
  
  if (!is.null(opt$motifs)) {
    message("Making chromavar analysis...")
    if (is.null(opt$cores)) {
      opt$cores <- 12
    }
    
    if (is.null(opt$assay)) {
      opt$assay <- "peaks" 
    }
    
    
    library(BiocParallel)
    register(MulticoreParam(opt$cores))
    
    if (startsWith(x = opt$motifs, prefix = "JASPAR")) {
      ## Chromvar
      #if (!require("JASPAR2024", quietly = TRUE)) BiocManager::install("JASPAR2024")
      
      library(JASPAR2022)
      
      ## Chromvar vertebrates
      # add motif information
      pfm <- getMatrixSet(
        x = JASPAR2022,
        opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
      )
      
      # add motif information
      combined.metacells <- AddMotifs(
        object = combined.metacells,
        genome = BSgenome.Mmusculus.UCSC.mm10,
        pfm = pfm,
        assay = opt$assay
      )
      DefaultAssay(combined.metacells) <- opt$assay
      
      combined.metacells <- RunChromVAR(
        object = combined.metacells,
        genome = BSgenome.Mmusculus.UCSC.mm10,
        new.assay.name = "chromvar_jaspar_vertebrates"
      )
      ## Chromvar mouse only
      
      # add new motif information
      pfm <- getMatrixSet(
        x = JASPAR2022,
        opts = list(collection = "CORE", species = 10090, all_versions = FALSE)
      )
      
      # add motif information
      combined.metacells <- AddMotifs(
        object = combined.metacells,
        genome = BSgenome.Mmusculus.UCSC.mm10,
        pfm = pfm,
        assay = opt$assay
      )
      DefaultAssay(combined.metacells) <- opt$assay
      
      combined.metacells <- RunChromVAR(
        object = combined.metacells,
        genome = BSgenome.Mmusculus.UCSC.mm10,
        new.assay.name = "chromvar_jaspar_mouse"
      )
    } else {
      if (endsWith(opt$motifs,suffix = ".rds"))
        pfm <- readRDS(opt$motifs)
      
      combined.metacells <- AddMotifs(
        object = combined.metacells,
        genome = BSgenome.Mmusculus.UCSC.mm10,
        pfm = pfm,
        assay = opt$assay
      )
      
      DefaultAssay(combined.metacells) <- opt$assay
      
      combined.metacells <- RunChromVAR(
        object = combined.metacells,
        genome = BSgenome.Mmusculus.UCSC.mm10,
        new.assay.name = "chromvar_pfm"
      )
      
    }
    

  }
  
  saveRDS(combined.metacells,paste0(opt$outdir,"/combined.metacells.with.macs.peak.rds"))
  
} else {
  saveRDS(peaks,paste0(opt$outdir,"/macs2.peaks.rds"))
  
}