library(Seurat)
library(dplyr)
library(getopt)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "sampleName",  "s", 1, "character", "sample name",
  "RNAcomp", "p", 1, "character", "range of components to consider for wnn analysis (eg 1:40 for RNA pca)",
  "ADTcomp", "q", 1, "character", "range of components to consider for wnn analysis (eg 1:50 for ADT pca)",
  "RNAnormalization", "r", 1, "character", "normalisation method for RNA (LogNormalize or SCT)",
  "SaveH5adFiles", "f", 1, "logical", "Save H5ad files for the different modalities"
  
), byrow=TRUE, ncol=5)

opt = getopt(spec)




if(is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "LogNormalize"
}

if(is.null(opt$SaveH5adFiles)) {
  opt$SaveH5adFiles <- FALSE
  
}


if (is.null(opt$RNAcomp)) {
  opt$RNAcomp <- c(1:40)
} else {
  ci <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][2])
  opt$RNAcomp <- c(ci:cf)
}

if (is.null(opt$ADTcomp)) {
  opt$ADTcomp <- c(1:50)
} else {
  ci <- as.numeric(strsplit(opt$ADTcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$ADTcomp,split = ":")[[1]][2])
  opt$ADTcomp <- c(ci:cf)
}

if (is.null(opt$outdir)) {
  opt$outdir <- "./"
}

print(opt)

dir.create(opt$outdir,recursive = T,showWarnings = F)

options(future.globals.maxSize = 120 * 1024 ^ 3) # for 120 Gb RAM


pbmc <- readRDS(opt$inputSeurat)
pbmc$orig.ident <- paste0(pbmc$donor,"_",pbmc$time)
pbmc <- pbmc[,pbmc$celltype.l2 != "Doublet"]
#pbmc <- pbmc[,pbmc$time == 0]  # first test on reference samples
pbmc$celltype.l1.5 <- pbmc$celltype.l1
pbmc$celltype.l1.5[pbmc$celltype.l1 == "other"] <- pbmc$celltype.l2[pbmc$celltype.l1 == "other"]

pbmc$celltype.l1.5[pbmc$celltype.l1 == "other T"] <- pbmc$celltype.l2[pbmc$celltype.l1 == "other T"]
gc()

pbmc.list <- SplitObject(pbmc, split.by = "orig.ident")

# integration of ADT data

pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
  DefaultAssay(x) <- "ADT"
  VariableFeatures(x) <- rownames(x[["ADT"]])
  x <- NormalizeData(x, normalization.method = 'CLR', margin = 2) %>% ScaleData(features = rownames(x)) %>% RunPCA(reduction.name = 'pca')
  return(x)})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = pbmc.list)
reference <- which(endsWith(names(pbmc.list),suffix = "_0"))

print(reference)

pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, 
                                       anchor.features = features, 
                                       reference = reference,
                                       reduction = "rpca",
                                       dims = opt$ADTcomp)

pbmc.combined <- IntegrateData(anchorset = pbmc.anchors,new.assay.name = "integratedADT")

Assays(pbmc.combined)


## Integration of SCT data

if (opt$RNAnormalization != "SCT") {
  rnaAssay = "RNA"
  pbmc.list <- lapply(X = pbmc.list, FUN = function(x) { DefaultAssay(x) <- "RNA";NormalizeData(x,normalization.method = "LogNormalize",assay = "RNA") %>% FindVariableFeatures(); return(x)})
  features <- SelectIntegrationFeatures(object.list = pbmc.list)
  pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
    x <- ScaleData(x,features = features);
    x <- RunPCA(x, features = features, verbose = FALSE);
    return(x)
  })
} else {
  rnaAssay = "SCT"
  pbmc.list <- lapply(X = pbmc.list, FUN = function(x) { DefaultAssay(x) <- "RNA";
  x <- SCTransform(x,conserve.memory = TRUE, verbose = FALSE); return(x)})
  features <- SelectIntegrationFeatures(object.list = pbmc.list)
  pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = features)
  
}

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features


pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, 
                                       anchor.features = features, 
                                       reduction = "rpca",
                                       reference = reference,
                                       normalization.method = opt$RNAnormalization,
                                       dims = opt$RNAcomp)

pbmc.combined.2 <- IntegrateData(anchorset = pbmc.anchors,normalization.method = opt$RNAnormalization,new.assay.name = paste0("integrated",rnaAssay))


pbmc.combined[["SCT"]] <- pbmc.combined.2[["SCT"]]

pbmc.combined[["integratedSCT"]] <- pbmc.combined.2[["integratedSCT"]]

remove(pbmc.combined.2)
gc()

Assays(pbmc.combined)
#saveRDS(pbmc.combined, paste0(opt$outdir,"/seuratCombined.rds"))

# WNN integration on integrated assay

DefaultAssay(pbmc.combined) <- 'integratedADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(pbmc.combined) <- rownames(pbmc.combined[["integratedADT"]]) 
pbmc.combined <- ScaleData(pbmc.combined)
pbmc.combined <- RunPCA(pbmc.combined,reduction.name = "apca")


if (opt$RNAnormalization == "SCT") {
  DefaultAssay(pbmc.combined) <- paste0("integrated",rnaAssay)
  pbmc.combined <- RunPCA(pbmc.combined) # integrated SCT slot can be used as scaled normalized counts
} else {
  DefaultAssay(pbmc.combined) <- paste0("integrated",rnaAssay)
  pbmc.combined <- ScaleData(pbmc.combined,features = features) %>% RunPCA() # integrated LogNormalized slot can be used as normalized counts
}

## Save in h5ad for SEACells
if (opt$SaveH5adFiles) {
SeuratDisk::SaveH5Seurat(pbmc.combined, filename =  paste0(opt$outdir,"/seurat.h5Seurat"),overwrite = T)
SeuratDisk::Convert(paste0(opt$outdir,"/seurat.h5Seurat"), dest =paste0(opt$outdir,"/seurat.ADT.h5ad"),assay ="integratedADT",overwrite = T)
SeuratDisk::Convert(paste0(opt$outdir,"/seurat.h5Seurat"), dest =paste0(opt$outdir,"/seurat.",rnaAssay,".h5ad"),assay =paste0("integrated",rnaAssay),overwrite = T)
system(command = paste0("rm -f ",paste0(opt$outdir,"/seurat.h5Seurat")))
}

pbmc.combined <- FindMultiModalNeighbors(
  pbmc.combined, 
  reduction.list = list("pca", "apca"), 
  dims.list = list(opt$RNAcomp, opt$ADTcomp), 
  modality.weight.name = "RNA.weight",
  return.intermediate = T
)
pbmc.combined <- RunUMAP(pbmc.combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")


pbmc.combined <- RunUMAP(pbmc.combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc.combined <- FindClusters(pbmc.combined, resolution = c(c(5:15)/10), graph.name = "wsnn", algorithm = 3)


saveRDS(pbmc.combined, paste0(opt$outdir,"/seuratCombinedWNN.rds"))

