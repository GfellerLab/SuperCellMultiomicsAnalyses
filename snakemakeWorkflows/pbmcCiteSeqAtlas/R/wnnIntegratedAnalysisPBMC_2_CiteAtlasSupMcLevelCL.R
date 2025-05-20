library(Seurat)
library(dplyr)
library(getopt)
library(SuperCellMultiomics)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'gamma', 'g', 1, "numeric", 'gamma to used for metacell identification',
  "sampleName",  "s", 1, "character", "sample name",
  "RNAcomp", "p", 1, "character", "range of components to consider for wnn analysis (eg 1:40 for RNA pca)",
  "ADTcomp", "q", 1, "character", "range of components to consider for wnn analysis (eg 1:50 for ADT pca)",
  "RNAnormalization", "r", 1, "character", "RNA normalisation method for integration  (only 'LogNormalize' at the moment)",
  "RNAmetacells", "m", 1, "character", "RNA normalisation method for metacells identification (default 'LogNormalize' or 'SCT' for SCTransform)",
  "SaveH5adFiles", "f", 1, "logical", "Save H5ad files for the different modalities"
  
), byrow=TRUE, ncol=5)

opt = getopt(spec)




if(is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "LogNormalize"
}

if(is.null(opt$RNAmetacells)) {
  opt$RNAmetacells <- "LogNormalize"
}

if(is.null(opt$SaveH5adFiles)) {
  opt$SaveH5adFiles <- FALSE
  
}


if (is.null(opt$RNAcomp)) {
  opt$RNAcomp <- c(1:30)
} else {
  ci <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][2])
  opt$RNAcomp <- c(ci:cf)
}

if (is.null(opt$ADTcomp)) {
  opt$ADTcomp <- c(1:30)
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


pbmc <- readRDS(opt$inputSeurat)
pbmc$orig.ident <- paste0(pbmc$donor,"_",pbmc$time)
pbmc <- pbmc[,pbmc$celltype.l2 != "Doublet"]
#pbmc <- pbmc[,pbmc$time == 0]  # first test on reference samples
pbmc$celltype.l1.5 <- pbmc$celltype.l1
pbmc$celltype.l1.5[pbmc$celltype.l1 == "other"] <- pbmc$celltype.l2[pbmc$celltype.l1 == "other"]

pbmc$celltype.l1.5[pbmc$celltype.l1 == "other T"] <- pbmc$celltype.l2[pbmc$celltype.l1 == "other T"]

gc()


#Per sample metacell identification
# 
# sampleNames <- c("P1_0","P2_0","P3_0","P4_0","P5_0","P6_0","P7_0","P8_0",
#                   "P1_2","P2_2","P3_2","P4_2","P5_2","P6_2","P7_2","P8_2",
#                   "P1_7","P2_7","P3_7","P4_7","P5_7","P6_7","P7_7","P8_7")
# 
# samplePaths = paste0("output/correlationAnalyzis/CITEseq/",sampleNames,"/singlecells_analysis/seuratWNN.rds")

pbmc$orig.ident <- paste0(pbmc$donor,"_",pbmc$time)

sampleNames <- unique(pbmc$orig.ident)

## test loop without a function!!

pbmc$celltype.l1.5 <- pbmc$celltype.l1
pbmc$celltype.l1.5[pbmc$celltype.l1 == "other"] <- pbmc$celltype.l2[pbmc$celltype.l1 == "other"]

pbmc$celltype.l1.5[pbmc$celltype.l1 == "other T"] <- pbmc$celltype.l2[pbmc$celltype.l1 == "other T"]

start_time <- Sys.time()

pbmc$orig.ident <- paste0(pbmc$donor,"_",pbmc$time)

sampleNames <- unique(pbmc$orig.ident)

## test loop without a function!!

pbmcMC_supL1.5.list <- list()

metaAll <- data.frame()

for (sampleName in sampleNames) {
  print(sampleName)
  print(dim(pbmc))
  pbmc.smp <- pbmc[,pbmc$orig.ident == sampleName & pbmc$celltype.l2 != "Doublet"]
  DefaultAssay(pbmc.smp) <- 'ADT'
  # we will use all ADT features for dimensional reduction
  # we set a dimensional reduction name to avoid overwriting the 
  VariableFeatures(pbmc.smp) <- rownames(pbmc.smp[["ADT"]])
  pbmc.smp <- NormalizeData(pbmc.smp, normalization.method = 'CLR', margin = 2) %>% 
    ScaleData() %>% RunPCA(reduction.name = 'apca')
  
  DefaultAssay(pbmc.smp) <- "RNA"
  
  if (opt$RNAmetacells == "LogNormalize") {
    pbmc.smp <- NormalizeData(pbmc.smp,normalization.method = "LogNormalize",assay = "RNA") %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
    
    pbmcMC <- SCimplify_for_Seurat(seurat = pbmc.smp,
                                   assay = c('RNA','ADT'),
                                   reduction = list("pca", "apca"), 
                                   dims = list(opt$RNAcomp, opt$ADTcomp),
                                   graph.name = "knn",
                                   label = "celltype.l1.5",
                                   kernel = T,
                                   gamma = opt$gamma)
    
  } else {
    pbmc.smp <-  SCTransform(pbmc.smp, vst.flavor = "v2", verbose = FALSE) %>%
      RunPCA(verbose = FALSE)
    
    pbmcMC <- SCimplify_for_Seurat(seurat = pbmc.smp,
                                   assay = c('SCT','ADT'),
                                   reduction = list("pca", "apca"), 
                                   dims = list(opt$RNAcomp, opt$ADTcomp),
                                   graph.name = "knn",
                                   label = "celltype.l1.5",
                                   kernel = T,
                                   gamma = opt$gamma)
  }
  
  
  remove(pbmc.smp)
  if (length(unique(pbmc$orig.ident)) > 1) {
    pbmc <- pbmc[,pbmc$orig.ident != sampleName]
  }
  gc() 
  
  pbmcMC_supL1.5.list[[sampleName]] <- pbmcMC
  
  
  
}

#saveRDS(pbmcMC_supL1.5.list,"pbmcMC_supL1.5.list.rds")


#pbmcMC_supL1.5.list <- readRDS("pbmcMC_supL1.5.list.rds") #identify with prior SCT normalization on RNA seq
library(STACAS)

# integration of ADT data

pbmcMC_supL1.5.list <- lapply(X = pbmcMC_supL1.5.list, FUN = function(x) {
  DefaultAssay(x) <- "ADT"
  VariableFeatures(x) <- rownames(x[["ADT"]])
  x <- RenameCells(x,add.cell.id = unique(x$orig.ident))
  x@misc$new.membership <- paste0(unique(x$orig.ident),"_", x@misc$membership)
  names(x@misc$new.membership) <- names(x@misc$membership)
  x <- NormalizeData(x, normalization.method = 'CLR', margin = 2) %>% ScaleData(features = rownames(x)) %>% RunPCA(reduction.name = 'pca')
  return(x)})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = pbmcMC_supL1.5.list)

reference <- which(endsWith(names(pbmcMC_supL1.5.list),suffix = "_0"))

print(reference)

pbmc.anchors <- FindAnchors.STACAS(object.list = pbmcMC_supL1.5.list, 
                                   anchor.features = features,
                                   cell.labels = "celltype.l1.5",
                                   reference = reference,
                                   dims = opt$ADTcomp)

pbmc.combined <- IntegrateData.STACAS(anchorset = pbmc.anchors,new.assay.name = "integratedADT")

pbmc.combined@misc$membership <- unlist(lapply(pbmcMC_supL1.5.list,FUN = function(x) {x@misc$new.membership}))

Assays(pbmc.combined)


## Integration of RNA data


rnaAssay = "RNA"
pbmcMC_supL1.5.list <- lapply(X = pbmcMC_supL1.5.list, FUN = function(x) { DefaultAssay(x) <- "RNA";x <- NormalizeData(x, verbose = FALSE); return(x)})
features <- SelectIntegrationFeatures(object.list = pbmcMC_supL1.5.list)
pbmcMC_supL1.5.list <- lapply(X = pbmcMC_supL1.5.list, FUN = function(x) {
  x <- ScaleData(x)
  x <- RunPCA(x, features = features, verbose = FALSE)
})




# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features



pbmc.anchors <- FindAnchors.STACAS(object.list = pbmcMC_supL1.5.list, 
                                   anchor.features = features,
                                   cell.labels = "celltype.l1.5",
                                   reference = reference,
                                   dims = opt$RNAcomp)

pbmc.combined.2 <- IntegrateData.STACAS(anchorset = pbmc.anchors,new.assay.name = paste0("integrated",rnaAssay))


pbmc.combined[[rnaAssay]] <- pbmc.combined.2[[rnaAssay]]

pbmc.combined[[paste0("integrated",rnaAssay)]] <- pbmc.combined.2[[paste0("integrated",rnaAssay)]]

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


DefaultAssay(pbmc.combined) <- paste0("integrated",rnaAssay)
pbmc.combined <- ScaleData(pbmc.combined)
pbmc.combined <- RunPCA(pbmc.combined,reduction.name = "pca")

pbmc.combined <- RunUMAP(pbmc.combined, dims = c(opt$RNAcomp), reduction = "pca", reduction.name = "rna.umap")
DimPlot(pbmc.combined,reduction = "rna.umap",group.by = "celltype.l1.5")

pbmc.combined <- RunUMAP(pbmc.combined, dims = c(opt$ADTcomp), reduction = "apca", reduction.name = "adt.umap")
DimPlot(pbmc.combined,reduction = "adt.umap",group.by = "celltype.l1.5")


pbmc.combined <- FindMultiModalNeighbors(
  pbmc.combined, 
  reduction.list = list("pca", "apca"), 
  dims.list = list(opt$RNAcomp, opt$ADTcomp), 
  modality.weight.name = "RNA.weight",
  return.intermediate = T
)
pbmc.combined <- RunUMAP(pbmc.combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

saveRDS(pbmc.combined, paste0(opt$outdir,"/seuratCombinedWNN_2.rds"))




















