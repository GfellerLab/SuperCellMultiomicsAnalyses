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
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "RNAnormalization", "n", 1, "character", "normalization method ('LogNormalization' or 'SCT')",
  "RNAcomp", "p", 1, "character", "range of components to consider for wnn analysis (eg 1:40 for RNA pca)",
  "ATACcomp", "q", 1, "character", "range of components to consider for wnn analysis (eg 1:50 for ATAC lsi)",
  "filterResolution", "f", 1, "numeric", "filtering resolution (default 0.5)."
), byrow=TRUE, ncol=5)

opt = getopt(spec)

if (is.null(opt$RNAcomp)) {
  opt$RNAcomp <- c(1:50)
} else {
  ci <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][2])
  opt$RNAcomp <- c(ci:cf)
}

if (is.null(opt$ATACcomp)) {
  opt$ATACcomp <- c(2:50)
} else {
  ci <- as.numeric(strsplit(opt$ATACcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$ATACcomp,split = ":")[[1]][2])
  opt$ATACcomp <- c(ci:cf)
}

if (is.null(opt$RNAnormalization)) {opt$RNAnormalization = "LogNormalize"}

if (is.null(opt$filterResolution)) {
  opt$filterResolution <- 0.5
} 


if (is.null(opt$outdir)) {opt$outdir = "output/10xMultiomeProstateHan22/"}

# Get the files
metacells.files <- list.files(opt$outdir,recursive = T,pattern = "metacells.rds",full.names = T,include.dirs = T)

metacell.objs <- lapply(metacells.files, FUN = function(x){ 
  x <- readRDS(x)
  x <- RenameCells(x,add.cell.id = unique(x$orig.ident))
  DefaultAssay(x) <- "RNA"
  x[["ATAC"]] <- NULL
  return(x)
}
)


#Same strategy as in the original study
#1 - integrate RNA,
#2 - discard low qual cluster
#3 - integrate ATAC/RNA on remaining cells
#4 - multimodal wnn

#Integrating RNA assay
metacell.objs <- lapply(X = metacell.objs, FUN = function(x) {
  x <- FindVariableFeatures(x)
})


## Integration of RNA data

if (opt$RNAnormalization != "SCT") {
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) { DefaultAssay(x) <- "RNA";NormalizeData(x,normalization.method = "LogNormalize",assay = "RNA") %>% FindVariableFeatures(); return(x)})
  features <- SelectIntegrationFeatures(object.list = metacell.objs)
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) {
    x <- ScaleData(x,features = features);
    x <- RunPCA(x, features = features, verbose = FALSE);
    return(x)
  })
  
  print("Finding anchors")
  
  anchors <- FindIntegrationAnchors(object.list = metacell.objs, 
                                    reduction = "rpca",
                                    anchor.features = features,
                                    dims = opt$RNAcomp)
  
  print("Integrating data")
  rna.All.combined <- IntegrateData(anchorset = anchors, dims = opt$RNAcomp)
  rm(metacell.objs,anchors)
  gc()
  
  DefaultAssay(rna.All.combined) <- "integrated"
  rna.All.combined <- ScaleData(rna.All.combined, verbose = FALSE)
  
} else {
  rnaAssay = "SCT"
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) { DefaultAssay(x) <- "RNA";x <- SCTransform(x, verbose = FALSE,vst.flavor = "v2"); return(x)})
  features <- SelectIntegrationFeatures(object.list = metacell.objs,nfeatures = 3000)
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) {
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  metacell.objs <- PrepSCTIntegration(object.list = metacell.objs, anchor.features = features)
  
  print("Finding anchors")
  
  anchors <- FindIntegrationAnchors(object.list = metacell.objs, 
                                    reduction = "rpca", 
                                    normalization.method = "SCT",
                                    anchor.features = features,
                                    dims = opt$RNAcomp)
  
  print("Integrating data")
  rna.All.combined <- IntegrateData(anchorset = anchors, 
                                    dims = opt$RNAcomp,
                                    normalization.method = "SCT")
  rm(metacell.objs,anchors)
  gc()
  
  DefaultAssay(rna.All.combined) <- "integrated"
}

saveRDS(rna.All.combined,"rna.all.combined.rds")

rna.All.combined <- RunPCA(rna.All.combined, verbose = FALSE)
rna.All.combined <- RunUMAP(rna.All.combined, dims = opt$RNAcomp)
rna.All.combined <- FindNeighbors(rna.All.combined, dims = opt$RNAcomp)
rna.All.combined <- FindClusters(rna.All.combined, resolution = c(1:10)*0.1)

genes <- c("Epcam","Krt8","Chga","Krt5","Svs5","Col1a2","Pecam1","Ptprc","S100a9","C1qa","Cd3e","Cd19","Plp1")
DefaultAssay(rna.All.combined) <- "SCT"
for (r in c(1:10)*0.1) {
  pdf(paste0(opt$outdir,'res',r,"_all_rna_mcs.pdf"))
  Idents(rna.All.combined) <- paste0('integrated_snn_res.',r)
  print(UMAPPlot(rna.All.combined,label = T))
  print(VlnPlot(rna.All.combined,features = c("nFeature_RNA"),pt.size = 0.001,log = T))
  print(VlnPlot(rna.All.combined,features = c("nCount_RNA"),pt.size = 0.001,log = T))
  print(FeaturePlot(rna.All.combined,features = genes[1:7]))
  print(FeaturePlot(rna.All.combined,features = genes[8:13]))
  
  
  dev.off()
}

DefaultAssay(rna.All.combined) <- "integrated"

filter.res <- paste0("integrated_snn_res.",opt$filterResolution)
# Remove low quality cluster
df <- rna.All.combined@meta.data %>%
  group_by_at(filter.res) %>%
  summarise(medNF = median(nFeature_RNA))

lowQualClust <- df[which.min(df$medNF),filter.res]

Idents(rna.All.combined) <- filter.res
rna.filtered.combined <- subset(rna.All.combined,idents = lowQualClust,invert = T)

remove(rna.All.combined)
gc()

# integrate RNA clean data
print("Integrate RNA clean data")
DefaultAssay(rna.filtered.combined) <- rnaAssay
metacell.objs <- SplitObject(rna.filtered.combined,split.by = "orig.ident")

if (opt$RNAnormalization != "SCT") {
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) { DefaultAssay(x) <- "RNA";NormalizeData(x,normalization.method = "LogNormalize",assay = "RNA") %>% FindVariableFeatures(); return(x)})
  features <- SelectIntegrationFeatures(object.list = metacell.objs)
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) {
    x <- ScaleData(x,features = features);
    x <- RunPCA(x, features = features, verbose = FALSE);
    return(x)
  })
  
  print("Finding anchors")
  
  anchors <- FindIntegrationAnchors(object.list = metacell.objs, 
                                    reduction = "rpca",
                                    anchor.features = features,
                                    dims = opt$RNAcomp)
  
  print("Integrating data")
  rna.filtered.combined <- IntegrateData(anchorset = anchors, dims = opt$RNAcomp)
  rm(metacell.objs,anchors)
  gc()
  
  DefaultAssay(rna.filtered.combined) <- "integrated"
  rna.filtered.combined <- ScaleData(rna.filtered.combined, verbose = FALSE)
  
} else {
  rnaAssay = "SCT"
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) { DefaultAssay(x) <- "RNA";x <- SCTransform(x, verbose = FALSE,vst.flavor = "v2"); return(x)})
  features <- SelectIntegrationFeatures(object.list = metacell.objs)
  metacell.objs <- lapply(X = metacell.objs, FUN = function(x) {
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  metacell.objs <- PrepSCTIntegration(object.list = metacell.objs, anchor.features = features)
  
  print("Finding anchors")
  
  anchors <- FindIntegrationAnchors(object.list = metacell.objs, 
                                    reduction = "rpca", 
                                    normalization.method = "SCT",
                                    anchor.features = features,
                                    dims = opt$RNAcomp)
  
  print("Integrating data")
  rna.filtered.combined <- IntegrateData(anchorset = anchors, 
                                         dims = opt$RNAcomp,
                                         normalization.method = "SCT")
  rm(metacell.objs,anchors)
  gc()
  
  DefaultAssay(rna.filtered.combined) <- "integrated"
}



rna.filtered.combined <- RunPCA(rna.filtered.combined, verbose = FALSE)
rna.filtered.combined <- RunUMAP(rna.filtered.combined, dims = opt$RNAcomp)
rna.filtered.combined <- FindNeighbors(rna.filtered.combined, dims = opt$RNAcomp)
rna.filtered.combined <- FindClusters(rna.filtered.combined, resolution = 0.2)


kept.cells <- as.vector(colnames(rna.filtered.combined))
print("Loading ATAC")
metacell.objs.atac <- lapply(metacells.files,function(x) {
  x <- readRDS(x)
  smp <- unique(x$orig.ident)
  print(smp)
  x <- RenameCells(x,add.cell.id = unique(x$orig.ident))
  x <- x[,colnames(x) %in% kept.cells]
  x[["RNA"]] <- NULL
  gc()
  return(x)
})

atac.filtered.combined <- merge(metacell.objs.atac[[1]],metacell.objs.atac[-1])

print("Creating combined object")
#################add ATAC assay to RNA assay
rna.filtered.combined[["ATAC"]]<- atac.filtered.combined[["ATAC"]]
rm(atac.filtered.combined)
gc()
DefaultAssay(rna.filtered.combined) <- "ATAC"
rna.filtered.combined <- RunTFIDF(rna.filtered.combined)
rna.filtered.combined <- FindTopFeatures(rna.filtered.combined, min.cutoff = 'q0')
rna.filtered.combined <- RunSVD(rna.filtered.combined)

DefaultAssay(rna.filtered.combined) <- "ATAC"
rna.filtered.combined <- harmony::RunHarmony(
  object = rna.filtered.combined,
  group.by.vars = 'orig.ident',
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE
)

#########re-compute the UMAP using corrected LSI embeddings
rna.filtered.combined <- RunUMAP(rna.filtered.combined, dims = opt$ATACcomp, reduction = 'harmony',reduction.name = "umap.harmony")
DefaultAssay(rna.filtered.combined) <- "ATAC"
rna.filtered.combined <- FindMultiModalNeighbors(rna.filtered.combined, reduction.list = list("pca", "harmony"), dims.list = list(opt$RNAcomp, opt$ATACcomp))
rna.filtered.combined <- FindClusters(rna.filtered.combined, resolution = 0.1,graph.name = "wsnn") 
rna.filtered.combined <- RunUMAP(rna.filtered.combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

p1 <- DimPlot(rna.filtered.combined, reduction = "umap",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(rna.filtered.combined, reduction = "umap.harmony",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(rna.filtered.combined,reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

pdf(file=paste0(opt$outdir,"/umaps.pdf"),width=1200,height=480)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()
## 

saveRDS(rna.filtered.combined,paste0(opt$outdir,"/combined.metacells.rds"))
write.csv(rna.filtered.combined@reductions[['wnn.umap']]@cell.embeddings,paste0(opt$outdir,"/wnn.umap.csv"))
write.csv(rna.filtered.combined@reductions[['umap.harmony']]@cell.embeddings,paste0(opt$outdir,"/wnn.umap.csv"))
write.csv(rna.filtered.combined@reductions[['wnn.umap']]@cell.embeddings,paste0(opt$outdir,"/wnn.umap.csv"))
write.csv(rna.filtered.combined@reductions[['wnn.umap']]@cell.embeddings,paste0(opt$outdir,"/wnn.umap.csv"))
write.csv(rna.filtered.combined@meta.data,paste0(opt$outdir,"/wnn.all.combined.meta.data.csv"))


















