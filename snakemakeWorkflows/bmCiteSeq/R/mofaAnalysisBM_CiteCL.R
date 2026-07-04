library(getopt)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds or name for laoding with SeuratData)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "python", "y", 1, "character", "python path with palantir installed to compute diff map",
  "n_factors", "n", 1, "integer", "number of factors for MOFA (edfault 30)"

), byrow=TRUE, ncol=5)

opt = getopt(spec)

library(reticulate)
use_python(opt$python)
library(Seurat)
library(dplyr)
library(MOFA2)




# opt <- list()
# opt$inputSeurat <- "input/bmCiteSeq/bmcite.rds"
# opt$outdir =  "output/bmCiteSeq/singlecells_analysis"
# opt$python =  "/opt/conda/envs/mofa_env/bin/python"    










# if help was asked, print a friendly message
# and exit with a non-zero error code
# test
#setwd("~/flamingo/beegfs01/tests/SuperCellMultiomicsAnalyses/")
#opt <- list()
#opt$inputSeurat <- "output/bmCiteSeq/singlecells_analysis/seuratWNN.rds"
#opt$outdir <- "output/bmCiteSeq/singlecells_analysis"
#opt$RNAcomp <- "1:30"
#opt$ATACcomp <- "1:18"
#opt$RNAnormalization <- "SCTransform"


if (is.null(opt$n_factors)) {
  opt$n_factors <- 30
}


if (is.null(opt$outdir)) {
  opt$outdir <- "./"
}

print(opt)

dir.create(opt$outdir,recursive = T,showWarnings = F)

if (endsWith(opt$inputSeurat,'.rds')) {
bm <- readRDS(opt$inputSeurat)
} else {
  bm <- SeuratData::LoadData(ds = opt$inputSeurat)

}



bm[["RNA"]] <- as(object = bm[["RNA"]], Class = "Assay")
bm[["ADT"]] <- as(object = bm[["ADT"]], Class = "Assay")

mofa <- create_mofa(bm, assays = c("RNA","ADT"))
mofa


model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- opt$n_factors 

mofa <- prepare_mofa(mofa,
                     model_options = model_opts
)

mofa <- run_mofa(mofa,outfile=  paste0(opt$outdir,"/mofa_model.hdf5"))


saveRDS(mofa, paste0(opt$outdir,"/MOFA.rds"))
r2 <- get_variance_explained(mofa)


# Variance explained per factor
r2_pf <- r2$r2_per_factor$group1

# Maximum variance explained across all views for each factor
max_r2 <- apply(r2_pf, 1, max)

# Keep factors explaining >= 1% variance in at least one view
threshold <- 0 # no thresholds
keep_factors <- names(max_r2[max_r2 >= threshold])

# Print retained factors
print(keep_factors)

Z <- get_factors(mofa, factors = keep_factors, groups = "all")
Z <- do.call(rbind, Z)

bm[["mofa"]] <- CreateDimReducObject(embeddings = Z,key = "MOFA")

adata <- anndata::AnnData(X = Matrix::t(GetAssayData(object = bm,
                                                     layer = "counts",assay = "RNA")),
			                             obs = bm@meta.data,
						     #raw = adata.raw,
						     obsm = list("X_mofa" = bm[["mofa"]]@cell.embeddings))

anndata::write_h5ad(adata,paste0(opt$outdir,"/adata_mofa.h5ad"))

#saveRDS(bm, paste0(opt$outdir,"/seuratMOFA.rds"))

