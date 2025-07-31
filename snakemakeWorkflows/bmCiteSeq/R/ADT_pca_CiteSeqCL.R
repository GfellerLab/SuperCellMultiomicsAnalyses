library(reticulate)
use_python("/opt/conda/envs/MetacellAnalysisToolkit/bin/python", required = TRUE)
library(Seurat)
library(dplyr)
library(getopt)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds or name for laoding with SeuratData)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "RNAnormalization", "r", 1, "character", "normalisation method for RNA (logNormalize or SCTransform)"

), byrow=TRUE, ncol=5)

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code
# test
# opt <- list()
# opt$fragmentFile <- "~/Documents/multiomicsMetacells/multiome_bm_data/fragments_files/bm_granulocyte_sorted_10k_atac_fragments.tsv.gz"
# opt$inputSeurat <- "bmMultiome"
# opt$outdir <- "output/correlationAnalyzis/bmMultiome/singlecell_analysis"
# frag.file <- opt$fragmentFile
# opt$RNAcomp <- "1:50"
# opt$ATACcomp <- "2:50"
# opt$minCutOff <- "q0"
# opt$RNAnormalization <- "SCTransform"



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

DefaultAssay(bm) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca')

adata <- anndata::AnnData(X = Matrix::t(GetAssayData(object = bm,slot = "counts",assay = "ADT")),
                          obs = bm@meta.data,
                          #raw = adata.raw,
                          obsm = list("X_apca" = bm[["apca"]]@cell.embeddings))

anndata::write_h5ad(adata,paste0(opt$outdir,"/seurat.ADT.h5ad"))
## Save in h5ad for SEACells
# SeuratDisk::SaveH5Seurat(bm, filename =  paste0(opt$outdir,"/seurat.h5Seurat"),overwrite = T)
# SeuratDisk::Convert(paste0(opt$outdir,"/seurat.h5Seurat"), dest =paste0(opt$outdir,"/seurat.ADT.h5ad"),assay ="ADT",overwrite = T)
# system(command = paste0("rm -f ",paste0(opt$outdir,"/seurat.h5Seurat")))
