library("igraph")
library("RANN")
library("WeightedCluster")
library("corpcor")
library("weights")
library("Hmisc")
library("Matrix")
library("patchwork")
library("plyr")
library("irlba")
library(Seurat)
library(SuperCellMultiomics)
library(SingleCellExperiment)
library(getopt)
#library(future.apply)
library(ggplot2)
#library(chromVAR)
# library(JASPAR2022)
library(TFBSTools)
#library(motifmatchr)
#library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(EnsDb.Hsapiens.v86)
# .libPaths(c("./input/", .libPaths()))
# options(timeout = 1000)
# SeuratData::InstallData("pbmcMultiome",lib = "./input/")
library(SeuratData)

## Loading ATAC data

spec = matrix(c(
  'gamma', 'g', 1, "numeric", "gamma to apply",
  'seed', 's', 1, 'numeric', "seed to apply",
  'inputRDS', "i",1, 'character', 'single cell rds object preprocessed',
  'ouptut', 'o',1, 'character', 'output file path'), 
  byrow=TRUE, ncol=5)

opt = getopt(spec)

print(opt)

gamma <- opt$gamma
seed <- opt$seed

#gamma <- as.numeric(commandArgs(TRUE)[1])
#seed <- as.numeric(commandArgs(TRUE)[2])

#InstallData("pbmcMultiome.SeuratData")
# data("pbmc.atac")


# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
# grange.counts <- StringToGRanges(rownames(pbmc.atac))
# grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
# pbmc.atac <- pbmc.atac[as.vector(grange.use), ]
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "hg38"
# 
# pbmc.atac <- CreateChromatinAssay(
#   counts = pbmc.atac@assays$ATAC@counts,
#   genome = 'hg38',
#   #fragments = frag.file,
#   annotation = annotations
# )

## Loading ATAC data and combining RNA and ATAC data in one seurat object with two assays



# data("pbmc.rna")
# 
# pbmc.rna[["ATAC"]] <- pbmc.atac
# 
# remove(pbmc.atac)
# 
# pbmc <- pbmc.rna
# 
# remove(pbmc.rna)
# 
# 
# pbmc <- pbmc[,pbmc$seurat_annotations != "filtered"]
# pbmc$coarse.annotation <- pbmc$seurat_annotations
# 
# #pbmc$coarse.annotation[grepl(pattern = "B",x = pbmc$coarse.annotation)] <- "B"
# 
# pbmc$coarse.annotation[grepl(pattern = "CD8 TEM",x = pbmc$coarse.annotation)] <- "CD8 Mem"
# 
# pbmc$coarse.annotation[grepl(pattern = "CD4 TEM",x = pbmc$coarse.annotation)] <- "CD4 Mem"
# pbmc$coarse.annotation[grepl(pattern = "CD4 TCM",x = pbmc$coarse.annotation)] <- "CD4 Mem"
# 
# Idents(pbmc) <- "coarse.annotation"
# pbmc[["CellName"]] <- rownames(pbmc@meta.data)
# 
# library(dplyr)
# options(future.globals.maxSize= 1000*1024^2)
# 
# DefaultAssay(pbmc) <- "RNA"
# pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
# 
# 
# 
# # Signac scATAC-seq workflow
# 
# 
# # ATAC analysis
# # We exclude the first dimension as this is typically correlated with sequencing depth
# #grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
# DefaultAssay(pbmc) <- "ATAC"
# pbmc <- RunTFIDF(pbmc)
# pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
# pbmc <- RunSVD(pbmc)
# pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
# 
# 
# # UMAP results
# 
# p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "coarse.annotation", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")+ NoLegend()
# p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "coarse.annotation", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")+ NoLegend()
# 
# p1 + p2
# 
# 
# ## Multimodal analysis with Seurat
# 
# 
# pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
# pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
# pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
# 
# 
# p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "coarse.annotation", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
# p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
# 
# p3
# 
# 
# 
# 
# p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "seurat_annotations", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
# p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "seurat_annotations", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
# p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "seurat_annotations", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
# p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))



detect_small_pop_multi <- function(seurat_object, 
                                   method, 
                                   seed, 
                                   cell_line_to_subset, 
                                   gamma, 
                                   cells_to_subset, 
                                   cells_to_keep){
  detect <- FALSE
  n_size <- 1
  while (detect == FALSE) {
    so.sb <- subset(seurat_object, CellName %in% c(cells_to_keep,cells_to_subset[1:n_size]))
    print(dim(so.sb))
    
    if ((method == "SC2_Kernel")){
      MC.v2_k <- seurat.mc.multi <- SCimplify_for_Seurat_v5(
        so.sb, assay = c('RNA','ATAC'),
        reduction = list("pca", "lsi"), 
        dims = list(1:40, 2:40), 
        graph.name = "knn",
        gamma = gamma
      )
      #print(unique(MC.v2_k@meta.data$annotation[MC.v2_k@meta.data$annotation_purity > 0.8]))
      if (cell_line_to_subset %in% MC.v2_k@meta.data$coarse.annotation[MC.v2_k@meta.data$coarse.annotation_purity > 0.8] | ncol(so.sb) == ncol(seurat_object)){
        detect <- TRUE
        results <- tibble::tibble(method = method, cell_type = cell_line_to_subset, pop_size = n_size, gamma = gamma, seed = seed)
      }
      else {
        detect <- FALSE
        n_size <- n_size+5
      }
    }
    
    else if ((method == "SC2_no_Kernel")){
      MC.v2_no_k <- seurat.mc.multi <- SCimplify_for_Seurat_v5(
        so.sb, assay = c('RNA','ATAC'),
        reduction = list("pca", "lsi"), 
        dims = list(1:40, 2:40), 
        graph.name = "knn",
        gamma = gamma,
        kernel = FALSE
      )
      #print(unique(MC.v2_no_k@meta.data$annotation[MC.v2_no_k@meta.data$annotation_purity > 0.8]))
      if (cell_line_to_subset %in% MC.v2_no_k@meta.data$coarse.annotation[MC.v2_no_k@meta.data$coarse.annotation_purity > 0.8] | ncol(so.sb) == ncol(seurat_object)){
        detect <- TRUE
        results <- tibble::tibble(method = method, cell_type = cell_line_to_subset, pop_size = n_size, gamma = gamma, seed = seed)
      }
      else {
        detect <- FALSE
        n_size <- n_size+5
      }
    }
  }
  return(results)
}

seurat_object <- readRDS(opt$inputRDS)

Idents(pbmc) <- "coarse.annotation"
res.df2_multi <- tibble::tibble(method = NA, cell_type = NA, pop_size = NA, gamma = NA, seed = NA)

for (cl in levels(pbmc)){
    print(cl)
    cells_of_interest <- rownames(subset(pbmc@meta.data, coarse.annotation == cl))
    other_cells <- rownames(subset(pbmc@meta.data, coarse.annotation != cl))
    
      set.seed(seed)
      shuffled_cells <- sample(cells_of_interest)
      res.sc2_k <- detect_small_pop_multi(pbmc,method = "SC2_Kernel",seed = seed,cl,gamma = gamma,cells_to_subset = shuffled_cells,cells_to_keep = other_cells)
      print(res.sc2_k)
      res.df2_multi <- rbind(res.df2_multi,res.sc2_k)
      res.sc2_no_k <- detect_small_pop_multi(pbmc,method = "SC2_no_Kernel",seed = seed,cl,gamma = gamma,cells_to_subset = shuffled_cells,cells_to_keep = other_cells)
      print(res.sc2_no_k)
      res.df2_multi <- rbind(res.df2_multi,res.sc2_no_k)
    
    
  }

res.df2_multi <- na.omit(res.df2_multi)
props <- c()
nb_cells_per_cl <- table(pbmc@meta.data$coarse.annotation)
for (j in 1:nrow(res.df2_multi)){
  cl_of_interest <- res.df2_multi[j,]$cell_type
  untouched_cl <- names(nb_cells_per_cl)[names(nb_cells_per_cl) != cl_of_interest]
  nb_untouched_cl <- sum(nb_cells_per_cl[untouched_cl])
  total_nb <- nb_untouched_cl+res.df2_multi[j,]$pop_size
  prop.tmp <-  res.df2_multi[j,]$pop_size / total_nb
  props <- c(props,prop.tmp*100)
}

res.df2_multi$proportion <- props

write.csv(na.omit(res.df2_multi),paste0(opt$output,"/","seed_",seed,".csv"))

