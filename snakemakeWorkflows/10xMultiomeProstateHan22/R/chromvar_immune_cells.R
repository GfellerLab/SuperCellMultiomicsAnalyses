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
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)


options(future.globals.maxSize= 8000*1024^2)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'input',     'i',1, "character", 'input RDS',
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "cores",      'c', 1, "character", "number of cores",
  "assay",      "a", 1, "character", "which chromatine assay to use for chromvar analysis",
  "motifs",    "p", 1, "character", "path to pfm matrix (.rds) or JASPAR version"
), byrow=TRUE, ncol=5)

opt = getopt(spec)

print(opt)

combined.metacells <- readRDS(opt$input)

if (is.null(opt$cores)) {
  opt$cores <- 12
}


library(BiocParallel)
register(MulticoreParam(opt$cores))

if (startsWith(x = opt$motifs, prefix = "JASPAR")) {
  ## Chromvar
  if (!require("JASPAR2024", quietly = TRUE)) BiocManager::install("JASPAR2024")
  
  library(JASPAR2024)
  
  ## Chromvar vertebrates
  # add motif information
  pfm <- getMatrixSet(
    x = JASPAR2024,
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
  )
  
  # add motif information
  combined.metacells <- AddMotifs(
    object = combined.metacells,
    genome = BSgenome.Mmusculus.UCSC.mm10,
    pfm = pfm,
    assay = opt$assay
  )
  ## Chromvar mouse only
  
  # add new motif information
  pfm <- getMatrixSet(
    x = JASPAR2024,
    opts = list(collection = "CORE", species = 10090, all_versions = FALSE)
  )
  
  # add motif information
  combined.metacells <- AddMotifs(
    object = combined.metacells,
    genome = BSgenome.Mmusculus.UCSC.mm10,
    pfm = pfm,
    assay = opt$assay
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

}

DefaultAssay(combined.metacells) <- opt$assay

combined.metacells <- RunChromVAR(
  object = combined.metacells,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  new.assay.name = "chromvar"
)



saveRDS(combined.metacells,
        paste0(opt$outdir,"/immune.combined.metacells.with.motifs.rds"))





