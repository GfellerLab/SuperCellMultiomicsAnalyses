library(Seurat)
library(dplyr)
# library(parallel)
library(dplyr)
library(scIntegrationMetrics) 
library(getopt)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'input',  'i', 1, "character", "REQUIRED : seurat object (.rds)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "RNAnormalization", "r", 1, "character", "RNA normalisation method for singlecell integration (default 'SCT')",
  "RNAcomp", "p", 1, "character", "range of RNA components to consider for integration integration  (eg opt$RNAcomp for RNA pca)",
  "ADTcomp", "q", 1, "character", "range of ATAC components to consider for integration integration (eg 2:50 for ATAC lsi)"  ), byrow=TRUE, ncol=5)

#test
#  opt <- list()
# setwd("./")
# opt$input <- "output/pbmcCiteSeqAtlas/pbmcCiteSeqAtlas_single_cells/seuratCombinedWNN.rds"
# opt$outdir <- "../workflow_level4_results/core_atlas/metacells/"
# opt$RNAnormalization <- "SCT"

# dir.create("../workflow_level4_results/core_atlas/metacells/",recursive = T)

opt = getopt(spec)

if (is.null(opt$umapPlot)) {
  opt$umapPlot <- F
}


if (!is.null(opt$RNAcomp)) {
  ci <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][2])
  opt$RNAcomp <- c(ci:cf)
} else {
  opt$RNAcomp <- 1:40
}



if (!is.null(opt$ADTcomp)) {
  ci <- as.numeric(strsplit(opt$ADTcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$ADTcomp,split = ":")[[1]][2])
  opt$ADTcomp <- c(ci:cf)
} else {
  opt$ADTcomp <- 1:50
}


if (is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "SCT"
}

rna.assay <- "RNA" # to get the HVG from the integrated assay
if(opt$RNAnormalization == "SCT") {
  rna.assay <- "SCT"
} else {
  rna.assay <- "RNA"
}

print(opt)

pbmc <- readRDS(opt$input)

pbmc
options(future.globals.maxSize = 80 * 1024 ^ 3) # for 50 Gb RAM





# Make unintegrated analysis



if (opt$RNAnormalization == "RNA") {
  
  DefaultAssay(pbmc) <- "RNA"
  if (class(pbmc[["RNA"]]) == "Assay5") {
    pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])
  } 
  
  VariableFeatures(pbmc) <- VariableFeatures(pbmc[[paste0('integrated',rna.assay)]])
  pbmc <- NormalizeData(pbmc) %>% ScaleData() %>% 
    RunPCA(reduction.name = "unintegrated_pca") %>% 
    RunUMAP(reduction = "unintegrated_pca",dims = opt$RNAcomp,
            reduction.name = "unintegrated_rna.umap")
} else {
  DefaultAssay(pbmc) <- "SCT"
  VariableFeatures(pbmc) <- VariableFeatures(pbmc[[paste0('integrated',rna.assay)]])
  pbmc <- RunPCA(pbmc, verbose = F, reduction.name = "unintegrated_pca")  %>% 
    RunUMAP(reduction = "unintegrated_pca",dims = opt$RNAcomp,
            reduction.name = "unintegrated_rna.umap")
}

DefaultAssay(pbmc) <- "ADT"
VariableFeatures(pbmc) <- VariableFeatures(pbmc[[paste0('integratedADT')]])
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc,reduction.name = "unintegrated_apca") %>% 
  RunUMAP(reduction = "unintegrated_apca",dims = opt$ADTcomp,
          reduction.name = "unintegrated_adt.umap")


dim.reds <- list("unintegrated_pca" = opt$RNAcomp,
                 "unintegrated_apca" = opt$ADTcomp,
                 "pca" =  opt$RNAcomp,
                 "apca"= opt$ADTcomp
)

metrics <- c('CiLISI', 'celltype_ASW')

res.all.sc.all.rep <- data.frame()
for (rep in c(1:5)) {
  set.seed(rep)
  sampled.cells <- sample(Cells(pbmc),0.25*ncol(pbmc))
  pbmc.sub <- subset(pbmc,cells = sampled.cells)
  
  metrics.celltype.l1.5 <- sapply(names(dim.reds),function(x){
    
    pbmc.sub[["reduction"]] <- Seurat::CreateDimReducObject(embeddings = (Seurat::Embeddings(pbmc.sub[[x]])[,dim.reds[[x]]]),
                                                        key = "component",
                                                        assay =  pbmc.sub[[x]]@assay.used)
    
    metrics <- scIntegrationMetrics::getIntegrationMetrics(pbmc.sub, 
                                                           metrics = metrics,
                                                           meta.label = "celltype.l1.5",
                                                           meta.batch = "orig.ident",
                                                           method.reduction = "reduction")
    return(metrics)
  },USE.NAMES = T)
  
  print("label 1.5 done")
  
  metrics.celltype.l2 <- sapply(names(dim.reds),function(x){
    
    pbmc.sub[["reduction"]] <- Seurat::CreateDimReducObject(embeddings = (Seurat::Embeddings(pbmc.sub[[x]])[,dim.reds[[x]]]),
                                                        key = "component",
                                                        assay =  pbmc.sub[[x]]@assay.used)
    
    metrics <- scIntegrationMetrics::getIntegrationMetrics(pbmc.sub, 
                                                           metrics = metrics,
                                                           meta.label = "celltype.l2",
                                                           meta.batch = "orig.ident",
                                                           method.reduction = "reduction")
    return(metrics)
  },USE.NAMES = T)
  print("label 2 done")
  
  res.df.celltype.l1.5 <- data.frame(t(metrics.celltype.l1.5))
  res.df <- data.frame(t(metrics.celltype.l2))
  
  res.df.celltype.l1.5$modality <- "RNA"
  res.df.celltype.l1.5$label <- "celltype.l1.5"
  res.df.celltype.l1.5$method <- rownames(res.df.celltype.l1.5)
  res.df.celltype.l1.5$modality[res.df.celltype.l1.5$method %in% c("apca","unintegrated_apca")] <- "Protein"
  res.df.celltype.l1.5$CiLISI <- unlist(res.df.celltype.l1.5$CiLISI)
  res.df.celltype.l1.5$celltype_ASW <- unlist(res.df.celltype.l1.5$celltype_ASW)
  res.df.celltype.l1.5$norm_cLISI <- unlist(res.df.celltype.l1.5$norm_cLISI)
  
  res.df <- data.frame(t(metrics.celltype.l2))
  res.df$label <- "celltype.l2"
  res.df$modality <- "Protein"
  res.df$method <- rownames(res.df)
  res.df$modality[res.df$method %in% c("apca","unintegrated_apca")] <- "Protein"
  res.df$CiLISI <- unlist(res.df$CiLISI)
  res.df$celltype_ASW <- unlist(res.df$celltype_ASW)
  res.df$norm_cLISI <- unlist(res.df$norm_cLISI)
  
  
  res.all.sc <- rbind(res.df.celltype.l1.5,res.df)
  res.all.sc$rep <- rep
  rbind(res.all.sc.all.rep,res.all.sc)
}

saveRDS(res.all.sc,paste0(opt$outdir,"/bench_res.rds"))























