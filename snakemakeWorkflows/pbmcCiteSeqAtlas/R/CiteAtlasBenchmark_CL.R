library(Seurat)
library(dplyr)
library(parallel)
library(SuperCell)
library(scIntegrationMetrics)
library(getopt)


spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds)",
  'inputSingleCells',  'j', 1, "character", "REQUIRED : we need sc anno to benchmark metacells seurat object (.rds)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "RNAcomp", "p", 1, "character", "range of components to consider for wnn analysis (eg 1:40 for RNA pca)",
  "ADTcomp", "q", 1, "character", "range of components to consider for wnn analysis (eg 1:50 for ADT pca)",
  "RNAnormalization", "r", 1, "character", "RNA normalisation method for integration  (default 'LogNormalize' or 'SCT' for SCTransform)"
), byrow=TRUE, ncol=5)

# opt <- list()
#  setwd("../../../SuperCellMultiomicsAnalyses/")
#
# opt$inputSeurat <- "output/pbmcCiteSeqAtlas/supMetacells_SCT_supStacas_lognorm_mean//g20/seuratCombinedWNN.rds"
# opt$inputSingleCells <- "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
# opt$RNAnormalization <- "RNA"
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

options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM


dir.create(opt$outdir,recursive = T,showWarnings = F)

rna.assay <- "RNA" # to get the HVG from the integrated assay
if(opt$RNAnormalization == "SCT") {
  rna.assay <- "SCT"
} else {
  rna.assay <- "RNA"
}


pbmc <- readRDS(opt$inputSeurat)

pbmc

# Make unintegrated analysis

options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM


if (opt$RNAnormalization != "SCT") {

  DefaultAssay(pbmc) <- "RNA"
  if (class(pbmc[["RNA"]]) == "Assay5") {
    pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])
  }

  VariableFeatures(pbmc) <- VariableFeatures(pbmc[[paste0('integrated',rna.assay)]])
  pbmc <- ScaleData(pbmc) %>%
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

meta.data.sc.all <- readRDS(opt$inputSingleCells)
meta.data.sc.all <- meta.data.sc.all@meta.data
meta.data.sc.all$orig.ident <- paste0(meta.data.sc.all$donor,"_",meta.data.sc.all$time)
meta.data.sc.all <- meta.data.sc.all[meta.data.sc.all$celltype.l2 != "Doublet",]
#pbmc <- pbmc[,pbmc$time == 0]  # first test on reference samples
meta.data.sc.all$celltype.l1.5 <- meta.data.sc.all$celltype.l1
meta.data.sc.all$celltype.l1.5[meta.data.sc.all$celltype.l1 == "other"] <- meta.data.sc.all$celltype.l2[meta.data.sc.all$celltype.l1 == "other"]

meta.data.sc.all$celltype.l1.5[meta.data.sc.all$celltype.l1 == "other T"] <- meta.data.sc.all$celltype.l2[meta.data.sc.all$celltype.l1 == "other T"]

gc()


rownames(meta.data.sc.all) <- paste0(meta.data.sc.all$orig.ident,".",rownames(meta.data.sc.all) )
gc()
n.cores <- 4
cl <- makeCluster(n.cores)

metrics.celltype.l1.5 <-  parSapply(cl,c(names(dim.reds)),
                                    FUN = function(x,
                                                   pbmc,
                                                   dim.reds,
                                                   meta.data.sc.all){

                                      library(scIntegrationMetrics)
                                      library(Seurat)
                                      library(SuperCell)

                                      #temporary fix on supercell_silhouette
                                      # source("config/SuperCellMultiomics/R/MetacellASW.R")
                                      metrics <- SuperCell:::getIntegrationMetricsMetacells(sobj.mc = pbmc,
                                                                                bio.label.col = "celltype.l1.5",
                                                                                batch.label.col = "orig.ident",
                                                                                dims = dim.reds[[x]],
                                                                                meta.data.sc =meta.data.sc.all[names(pbmc@misc$membership),],
                                                                                reduction.name = x)
                                    },
                                    pbmc = pbmc,
                                    dim.reds =dim.reds,
                                    meta.data.sc.all = meta.data.sc.all,USE.NAMES = T)

stopCluster(cl)

n.cores <- 4
cl <- makeCluster(n.cores)

metrics.celltype.l2 <- parSapply(cl,c(names(dim.reds)),
                                 FUN = function(x,
                                                pbmc,
                                                dim.reds,
                                                meta.data.sc.all){

                                   library(scIntegrationMetrics)
                                   library(Seurat)
                                   library(SuperCell)
                                   # source("config/SuperCellMultiomics/R/MetacellASW.R")

                                   metrics <- SuperCell:::getIntegrationMetricsMetacells(sobj.mc = pbmc,
                                                                             bio.label.col = "celltype.l2",
                                                                             batch.label.col = "orig.ident",
                                                                             dims = dim.reds[[x]],
                                                                             meta.data.sc =meta.data.sc.all[names(pbmc@misc$membership),],
                                                                             reduction.name = x)
                                   return(metrics)
                                 },
                                 pbmc = pbmc,
                                 dim.reds =dim.reds,
                                 meta.data.sc.all = meta.data.sc.all,USE.NAMES = T)
stopCluster(cl)



metrics.celltype.l1.5 <- data.frame(t(metrics.celltype.l1.5))
res.df <- data.frame(t(metrics.celltype.l2))

metrics.celltype.l1.5$modality <- "RNA"
metrics.celltype.l1.5$label <- "celltype.l1.5"
metrics.celltype.l1.5$method <- rownames(metrics.celltype.l1.5)
metrics.celltype.l1.5$modality[metrics.celltype.l1.5$method %in% c("apca","unintegrated_apca")] <- "Protein"
metrics.celltype.l1.5$CiLISI <- unlist(metrics.celltype.l1.5$CiLISI)
metrics.celltype.l1.5$celltype_ASW <- unlist(metrics.celltype.l1.5$celltype_ASW)
metrics.celltype.l1.5$norm_cLISI <- unlist(metrics.celltype.l1.5$norm_cLISI)

res.df <- data.frame(t(metrics.celltype.l2))
res.df$label <- "celltype.l2"
res.df$modality <- "RNA"
res.df$method <- rownames(res.df)
res.df$modality[res.df$method %in% c("apca","unintegrated_apca")] <- "Protein"
res.df$CiLISI <- unlist(res.df$CiLISI)
res.df$celltype_ASW <- unlist(res.df$celltype_ASW)
res.df$norm_cLISI <- unlist(res.df$norm_cLISI)


res.all.sc <- rbind(metrics.celltype.l1.5,res.df)

saveRDS(res.all.sc,paste0(opt$outdir,"/bench_res.rds"))
