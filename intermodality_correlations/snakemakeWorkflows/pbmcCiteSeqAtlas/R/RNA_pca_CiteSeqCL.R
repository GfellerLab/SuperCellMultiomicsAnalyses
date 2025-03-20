library(Seurat)
library(dplyr)
library(getopt)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds or name for laoding with SeuratData)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  #"sampleName", "s", 1, "character", "sample name",
  "RNAnormalization", "r", 1, "character", "normalisation method for RNA (logNormalize or SCTransform)"
  
  
), byrow=TRUE, ncol=5)

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code
# test
# opt <- list()
# opt$fragmentFile <- "~/Documents/multiomicsMetacells/multiome_pbmc_data/fragments_files/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
# opt$inputSeurat <- "pbmcMultiome"
# opt$outdir <- "output/correlationAnalyzis/pbmcMultiome/singlecell_analysis"
# frag.file <- opt$fragmentFile 
# opt$RNAcomp <- "1:50"
# opt$ATACcomp <- "2:50"
# opt$minCutOff <- "q0"
# opt$RNAnormalization <- "SCTransform"

if(is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "logNormalize"
  
}



if (is.null(opt$outdir)) {
  opt$outdir <- "./"
}

print(opt)

dir.create(opt$outdir,recursive = T,showWarnings = F)

if (endsWith(opt$inputSeurat,'.rds')) {
  pbmc <- readRDS(opt$inputSeurat)
  pbmc[["weighted.nn"]] <- NULL
  pbmc[["wknn"]] <- NULL
  pbmc[["wsnn"]] <- NULL
  pbmc[["pca"]] <- NULL
  pbmc[["apca"]] <- NULL
  pbmc[["wnn.umap"]] <- NULL
  # pbmc$orig.ident <- paste0(pbmc$donor,"_",pbmc$time)
  #pbmc <- pbmc[,pbmc$orig.ident == opt$sampleName]
  #pbmc <- pbmc[,pbmc$celltype.l2 != "Doublet"]
  #pbmc <- pbmc[,pbmc$time == 0]  # first test on reference samples
  # pbmc$celltype.l1.5 <- pbmc$celltype.l1
  # pbmc$celltype.l1.5[pbmc$celltype.l1 == "other"] <- pbmc$celltype.l2[pbmc$celltype.l1 == "other"]
  # 
  # pbmc$celltype.l1.5[pbmc$celltype.l1 == "other T"] <- pbmc$celltype.l2[pbmc$celltype.l1 == "other T"]
  
  gc()
} else {
  pbmc <- SeuratData::LoadData(ds = opt$inputSeurat)

}


DefaultAssay(pbmc) <- "RNA"
pbmc[["ADT"]] <- NULL
gc()

if (opt$RNAnormalization == "SCT") {
  rnaAssay = "SCT"
  pbmc <- SCTransform(pbmc, verbose = FALSE,conserve.memory = T) %>% RunPCA() 
} else {
  rnaAssay = "RNA"
  pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",assay = "RNA") %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
}

adata <- anndata::AnnData(X = Matrix::t(GetAssayData(object = pbmc,slot = "counts",assay = "RNA")),
                          obs = pbmc@meta.data,
                          #raw = adata.raw,
                          obsm = list("X_pca" = pbmc[["pca"]]@cell.embeddings))

anndata::write_h5ad(adata,paste0(opt$outdir,"/seurat.RNA.h5ad"))

## Save in h5ad for SEACells
# SeuratDisk::SaveH5Seurat(pbmc, filename =  paste0(opt$outdir,"/seurat.h5Seurat"),overwrite = T)
# SeuratDisk::Convert(paste0(opt$outdir,"/seurat.h5Seurat"), dest =paste0(opt$outdir,"/seurat.",rnaAssay,".h5ad"),assay =rnaAssay,overwrite = T)
# if (rnaAssay != "RNA") {
# SeuratDisk::Convert(paste0(opt$outdir,"/seurat.h5Seurat"), dest =paste0(opt$outdir,"/seurat.RNA.h5ad"),assay ="RNA",overwrite = T)
# }
# system(command = paste0("rm -f ",paste0(opt$outdir,"/seurat.h5Seurat")))


