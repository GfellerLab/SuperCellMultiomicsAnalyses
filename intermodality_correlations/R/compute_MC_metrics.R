
# Libraries ---------------------------------------------------------------
library(Seurat)
library(Signac)
library(getopt)
library(future)
library(presto)
library(BiocParallel)
library(SuperCellMultiomics)


# Functions ---------------------------------------------------------------
source("R/functions/metacellCompactnessSeparation.R")

# Parameters --------------------------------------------------------------

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'mcSeurat',  'm', 1, "character", "Path to metacell data",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "scSeurat", "s", 1, "character", "Path to singlecell data",
  "pythonSeacellEnv", "e", 1, "character", "path to python",
  "RNAcomp", "n", 1, "character", "range of RNA components to consider for metacell identification  (eg 1:50 for RNA pca)",
  "ATACcomp", "q", 1, "character", "range of ATAC components to consider for metacell identification (eg 2:50 for ATAC pca)"
), byrow=TRUE, ncol=5)

opt = getopt(spec)

library(reticulate)
use_python(opt$pythonSeacellEnv)
library(MetacellAnalysisToolkit)


dir.create(opt$outdir, recursive = T, showWarnings = F)
print(opt)


# Load single-cell data  ------------------------------
seurat.sc <- readRDS(opt$scSeurat)

# Load seurat object with metacells ---------------------------------------

seurat.obj <- readRDS(opt$mcSeurat)


###########################################################################
###########################################################################
###                                                                     ###
###                     COMPUTE METACELL QC METRICS                     ###
###                                                                     ###
###########################################################################
###########################################################################

# compute metacell QC metrics ---------------------------------------------

if (!is.null(opt$RNAcomp)) {
  ci <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][2])
  opt$RNAcomp <- c(ci:cf)
}

if (!is.null(opt$ATACcomp)) {
  ci <- as.numeric(strsplit(opt$ATACcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$ATACcomp,split = ":")[[1]][2])
  opt$ATACcomp <- c(ci:cf)
}

membership_df <- data.frame("membership" = paste0("Metacell_", na.exclude(seurat.obj@misc$membership)),
                            row.names = colnames(seurat.sc)[!is.na(seurat.obj@misc$membership)])

DefaultAssay(seurat.sc) <- "RNA"
pca_diffusion_comp <- get_diffusion_comp(sc.obj = seurat.sc, dims = opt$RNAcomp)
colnames(pca_diffusion_comp) <- c(1:ncol(pca_diffusion_comp))
pca_diffusion_comp <-  as.matrix(pca_diffusion_comp)
seurat.sc[["pca_diffusion"]] <- CreateDimReducObject(embeddings = pca_diffusion_comp, assay = "RNA",key = "DM_")

lsi_diffusion_comp <- get_diffusion_comp(sc.obj = seurat.sc, dims = opt$ATACcomp, sc.reduction = "lsi")
colnames(lsi_diffusion_comp) <- c(1:ncol(lsi_diffusion_comp))
lsi_diffusion_comp <-  as.matrix(lsi_diffusion_comp)
seurat.sc[["lsi_diffusion"]] <- CreateDimReducObject(embeddings = lsi_diffusion_comp, assay = "ATAC", key = "ADM_")

seurat.obj$compactness_pca <- mc_compactness(cell.membership = membership_df,
                                             sc.obj = seurat.sc,
                                             sc.reduction = "pca_diffusion")
seurat.obj$compactness_lsi <- mc_compactness(cell.membership = membership_df,
                                             sc.obj = seurat.sc,
                                             sc.reduction = "lsi_diffusion")
seurat.obj$separation_pca <- mc_separation(cell.membership = membership_df,
                                           sc.obj = seurat.sc,
                                           sc.reduction = "pca_diffusion")
seurat.obj$separation_lsi <- mc_separation(cell.membership = membership_df,
                                           sc.obj = seurat.sc,
                                           sc.reduction = "lsi_diffusion")

message("computing ASW")
seurat.obj$asw_lsi_diffusion <- mc_ASW(seurat.mc = seurat.obj,
                                       seurat.sc = seurat.sc,
                                       cell.membership = membership_df,
                                       reduction.name.sc = "lsi_diffusion")

seurat.obj$asw_pca_diffusion <- mc_ASW(seurat.mc  = seurat.obj,
                                       seurat.sc = seurat.sc,
                                       cell.membership = membership_df,
                                       reduction.name.sc = "pca_diffusion")


message("all metrics computed")
saveRDS(seurat.obj, file = paste0(opt$outdir, "/seurat.multiome.mcMetrics.rds"))
