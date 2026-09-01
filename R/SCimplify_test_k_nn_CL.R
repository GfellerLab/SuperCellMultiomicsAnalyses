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

# pct.to.hide <- c(0:10)*10
k.to.test <- c(10,20,30,40,60,80,100)
meta.data.all <- data.frame()

for (k in k.to.test) {

  if(!is.null(opt$Mod2assay)) {
  seurat.mc.multi <- SCimplify_for_Seurat_v5(
    seurat,
    assay = c(opt$Mod1assay,opt$Mod2assay),
    reduction = list(opt$Mod1red,opt$Mod2red),
    dims = list(opt$Mod1comp, opt$Mod2comp),
    graph.name = "knn",
    k.knn = k,
    label = opt$guiding_label,
    gamma = opt$gamma
  )
  
  } else {
    seurat.mc.multi <- SCimplify_for_Seurat_v5(
      seurat,
      assay = c(opt$Mod1assay),
      reduction = list(opt$Mod1red),
      dims = list(opt$Mod1comp),
      k.knn=k,
      label = opt$guiding_label,
      gamma = opt$gamma
    )
  }

  saveRDS(seurat.mc.multi@misc,paste0(opt$outdir,"/k_",k,"_nn_memberships.rds"))

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
  meta.data$k_nn <- k
  meta.data.all <- rbind(meta.data.all,meta.data)

}

write.csv(meta.data.all,paste0(opt$outdir,"/results_test_k_nn.csv"))

