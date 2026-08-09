library(SeuratData)
library(Seurat)
library(Signac)
library(dplyr)
library(SuperCell)
library(getopt)
library(MetacellAnalysisToolkit)

# devtools::install_github("GfellerLab/SuperCellMultiomics",auth_token = "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")
# iterative_sampling <- function(sc_data, labels_col,percentages,seed=123) {
#   set.seed(seed)
#
#   initial_labels <- rownames(sc_data[[labels_col]])
#   total_sc <- length(initial_labels)
#   labels_to_hide <- list()
#   sampled_labels <- c()
#   # Iteratively sample based on percentages
#   for (percentage in percentages) {
#     # each time sample out labels from the ones retained in previous iteration
#     remaining_labels <- setdiff(initial_labels, sampled_labels)
#     sampled_labels <- c(sampled_labels,sample(remaining_labels, size = as.integer(percentage * total_sc / 100)))
#     labels_to_hide <- append(labels_to_hide, list(sampled_labels))
#   }
#
#   return(labels_to_hide)
# }


spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds) or dataset name (from SeuratData package)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "Mod1comp", "p", 1, "character", "range of Mod1 components to consider for metacell identification  (eg 1:50 for Mod1 pca)",
  "Mod2comp", "q", 1, "character", "range of Mod2 components to consider for metacell identification (eg 2:50 for Mod2 pca)",
  "Mod1red", "v", 1, "character", "Dim red name for assay 1",
  "Mod2red", "w", 1, "character", "Dim red name for assay 2",
  "Mod1assay", "r", 1, "character", "Mod1 assay name with computed pca (default Mod1)",
  "Mod2assay", "a", 1, "character", "Mod2 assay name with computed lsi (default Mod2)" ,
  "gamma", "g", 1, "numeric", "gamma used for metacell identificaiton",
  "seed", "s", 1, "numeric", "seed used for sampling labels to hide",
  "guiding_label", "l", 1, "character", "guiding label col"

), byrow=TRUE, ncol=5)

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code
# test
# opt <- list()
# opt$aggregateFragmentfile <-T
# opt$inputSeurat <- "../../output/bmCiteSeq/singlecells_analysis/seuratWNN.rds"
# # opt$outdir <- "output/correlationAnalyzis/pbmcMultiome/randomMetacells/"
# opt$Mod1comp <- "1:30"
# # opt$Mod2comp <- "1:18"
# # opt$aggregateFragmentfile <- TRUE
# opt$gamma <- 75
# opt$Mod1assay <- "RNA"
# opt$Mod1red <- "pca"
#
#
# opt$seed <- 2025
# opt$guiding_label <- "celltype.l2"


# if (is.null(opt$Mod1assay)) {
#   opt$Mod1assay <- "SCT"
# }
#
# if (is.null(opt$Mod2assay)) {
#   opt$Mod2assay <- "ATAC"
# }
#
#
# if (is.null(opt$Mod1red)) {
#   opt$Mod1red <- "pca"
# }
#
# if (is.null(opt$Mod2red)) {
#   opt$Mod2red <- "lsi"
# }


if (!is.null(opt$Mod1comp)) {
  ci <- as.numeric(strsplit(opt$Mod1comp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$Mod1comp,split = ":")[[1]][2])
  opt$Mod1comp <- c(ci:cf)
}

if (!is.null(opt$Mod2comp)) {
  ci <- as.numeric(strsplit(opt$Mod2comp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$Mod2comp,split = ":")[[1]][2])
  opt$Mod2comp <- c(ci:cf)
}

if (is.null(opt$outdir)) {
  opt$outdir <- "./"
}



print(opt)

dir.create(opt$outdir,recursive = T,showWarnings = F)


seurat <- readRDS(opt$inputSeurat)

seurat.pseudobulk <- AggregateExpression(seurat,group.by = "celltype.l2",slot = "counts",return.seurat = T)
seurat.pseudobulk <- NormalizeData(seurat.pseudobulk)

distances <- dist(as.matrix(Matrix::t(GetAssayData(seurat.pseudobulk,assay = "RNA"))))

library(reshape2)

distances.df <- melt(as.matrix(distances), varnames = c("row", "col"))


distances.df

getClosestSubtype <- function(distances) {
  res <- as.character(rownames(as.matrix(distances))[order(as.matrix(distances)[,1])[2]])
  return(res)
}



closestSubtypes <- apply(data.frame(as.matrix(distances)),MARGIN = 2,FUN = getClosestSubtype)
names(closestSubtypes) <- rownames(as.matrix(distances))

closestSubtypes <- data.frame(closestSubtypes)
rownames(closestSubtypes) <- gsub(x=rownames(closestSubtypes),pattern = "-",replacement = "_")
closestSubtypes$closestSubtypes <- gsub(x=closestSubtypes$closestSubtypes,pattern = "-",replacement = "_")

write.csv(closestSubtypes,paste0(opt$outdir,"/closest_celltypes.csv"))

# pct.to.hide <- c(0:10)*10
pct.to.hide <- c(0,5,10,25,50,NA)
meta.data.all <- data.frame()

for (p in pct.to.hide) {
  if (!is.na(p)) {
  guiding_label_col <- paste0("guiding_label_",p,"_pct_shuffled")
  
  set.seed(opt$seed)
  seurat[["label_col"]] <- seurat[[opt$guiding_label]]
  cells.to.shuffle <- sample(colnames(seurat),
                             size = as.integer(pct.to.shuffle[p] * ncol(seurat) / 100))
  seurat@meta.data[cells.to.shuffle,"label_col"] <- plyr::mapvalues(x=seurat@meta.data[cells.to.shuffle,"label_col"],
                                                                    from = rownames(closestSubtypes),
                                                                    to = closestSubtypes$closestSubtypes)
  
  } else {
    guiding_label_col <- "unsupervised"
    seurat[["label_col"]] <- "None"
  }

  

  if(!is.null(opt$Mod2assay)) {
  seurat.mc.multi <- SCimplify_for_Seurat_v5(
    seurat,
    assay = c(opt$Mod1assay,opt$Mod2assay),
    reduction = list(opt$Mod1red,opt$Mod2red),
    dims = list(opt$Mod1comp, opt$Mod2comp),
    graph.name = "knn",
    label = "label_col",
    gamma = opt$gamma
  )
  
  } else {
    seurat.mc.multi <- SCimplify_for_Seurat_v5(
      seurat,
      assay = c(opt$Mod1assay),
      reduction = list(opt$Mod1red),
      dims = list(opt$Mod1comp),
      label = "label_col",
      gamma = opt$gamma
    )
  }

  saveRDS(seurat.mc.multi,paste0(opt$outdir,"/",guiding_label_col,"seurat.mc.rds"))

  reds <- Reductions(seurat)[endsWith(x = Reductions(seurat),suffix = "_diffusion")]
  # we discard metacell outliers with na.exclude
  membership_df <- data.frame("membership" = paste0("Metacell_",na.exclude(seurat.mc.multi@misc$membership)),
                              row.names = colnames(seurat)[!is.na(seurat.mc.multi@misc$membership)])

  for (r in reds) {

    red.ori <- strsplit(r,"_")[[1]][1]

  seurat.mc.multi[[paste0("compactness_",red.ori)]] <- mc_compactness(cell.membership = membership_df,
                                                                           sc.obj = seurat,
                                                                           sc.reduction = r)
  seurat.mc.multi[[paste0("separation_",red.ori)]] <- mc_separation(cell.membership = membership_df,
                                                                       sc.obj = seurat,
                                                                       sc.reduction =r)


  }





  meta.data <- seurat.mc.multi@meta.data
  meta.data$guiding_labels <- guiding_label_col
  meta.data.all <- rbind(meta.data.all,meta.data)

}

write.csv(meta.data.all,paste0(opt$outdir,"/results_test_shuffle_closest.csv"))

# library(ggplot2)
# library(Cairo)
# options(bitmapType = "cairo")
# ggplot(meta.data.all,aes(x=as.factor(pct_known),y=celltype.l2_purity)) + geom_boxplot()
#
# ggplot(meta.data.all,aes(x=as.factor(pct_known),y=size)) + geom_boxplot()
#
# ggplot(meta.data.all,aes(x=as.factor(pct_known),y=compactness_pca)) + geom_boxplot() + scale_y_log10()
#
# ggplot(meta.data.all,aes(x=as.factor(pct_known),y=separaton_pca)) + geom_boxplot() + scale_y_log10()
#
# ggplot(meta.data.all,aes(x=as.factor(pct_known),y=compactness_apca)) + geom_boxplot() + scale_y_log10()
#
# ggplot(meta.data.all,aes(x=as.factor(pct_known),y=separaton_apca)) + geom_boxplot() + scale_y_log10() + ggsignif::geom_signif()
#
# DefaultAssay(seurat.mc.multi) <- "ADT"
#
# seurat.mc.multi <- NormalizeData(seurat.mc.multi,normalization.method = "CLR",margin = 2)
# FeatureScatter.SuperCell(seurat.mc.multi[,seurat.mc.multi$celltype.l1 %in% "T cell"],feature1 = "adt_CD8a",feature2 = "adt_CD4")
#
# DimPlotSC(seurat = seurat,
#           seurat.mc = seurat.mc.multi,
#           metacell.col = "celltype.l2",
#           sc.col = "celltype.l2") + theme_classic()
