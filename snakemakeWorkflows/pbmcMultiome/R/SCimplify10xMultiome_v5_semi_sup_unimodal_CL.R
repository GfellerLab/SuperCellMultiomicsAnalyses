library(SeuratData)
library(Seurat)
library(Signac)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SuperCellMultiomics)
library(getopt)
library(doParallel)


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
  "Mod1assay", "r", 1, "character", "Mod1 assay name with computed pca (default Mod1)",
  "Mod2assay", "a", 1, "character", "Mod2 assay name with computed lsi (default Mod2)" , 
  "k.wnn", "k", 1, "numeric", "k for the knn used in the wnn used for metacell identification",
  "gamma", "g", 1, "numeric", "gamma used for metacell identificaiton",
  "kernel", "e", 1, "logical", "whether to use a kernel",
  "randomMetacells", "d", 1, "logical","construct random metacell at the specified gamma",
  "memberships", "c", 1 , "character", "memberships already computed (e.g. with SEACells) csv table: 1st column cell name, 2nd column metacell name",
  "prefixMembership", "x", 1, "character", "metacell name prefix memebership to disccard (e.g. SEACell-)",
  "inputSeuratMetacell", "m", 1, "character", "seurat metacell object if available to rescale directly",
  "aggregateFragmentfile", "f", 0, "logical", 'whether to aggregate fragment file or not (default FALSE)',
  "returnMemberships", "b", 0, "logical", "wether to return only memberships (default seurat object with memberships inn misc slot)"
  
), byrow=TRUE, ncol=5)

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code
# test
opt <- list()
# opt$aggregateFragmentfile <-T
opt$inputSeurat <- "../../output/bmCiteSeq/singlecells_analysis/seuratWNN.rds/"
# opt$outdir <- "output/correlationAnalyzis/pbmcMultiome/randomMetacells/"
opt$Mod1comp <- "1:40"
opt$Mod2comp <- "2:40"
opt$kernel <- TRUE
# opt$aggregateFragmentfile <- TRUE
opt$gamma <- 20
opt$Mod1assay <- "SCT"
opt$seed <- 2025
opt$guiding_label <- "seurat_annotations"
pt$pythonSeacellEnv <- 
  
  
  if (is.null(opt$Mod1assay)) {
    opt$Mod1assay <- "RNA"
  }

if (is.null(opt$Mod2assay)) {
  opt$Mod2assay <- "ATAC"
}


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

if (is.null(opt$randomMetacells)) {
  opt$randomMetacells <- F
}

print(opt)

dir.create(opt$outdir,recursive = T,showWarnings = F)


seurat <- readRDS(opt$inputSeurat)
pct.to.hide <- c(0:10)*10
pct.to.hide <- c(0,25,50,75,100)
guiding_label_cols <- paste0("guiding_label_",pct.to.hide,"_pct_hidden")
names(pct.to.hide) <- guiding_label_cols
meta.data.all <- data.frame()

for (p in guiding_label_cols) {
  set.seed(opt$seed)
  seurat[["label_col"]] <- seurat[[opt$guiding_label]]
  
  if (p != "guiding_label_100_pct_hidden") {
    cells.to.hide <- sample(colnames(seurat),
                            size = as.integer(pct.to.hide[p] * ncol(seurat) / 100))
    seurat@meta.data[cells.to.hide,"label_col"] <- NA
  } 
  seurat.mc.multi <- SCimplify_for_Seurat_v5(
    seurat, 
    assay = c(opt$Mod1assay,opt$Mod2assay),
    reduction = list("pca", "lsi"), 
    dims = list(opt$Mod1comp, opt$Mod2comp), 
    graph.name = "knn",
    label = "label_col",
    gamma = opt$gamma
  )
  
  
  library(MetacellAnalysisToolkit)
  # we discard metacell outliers with na.exclude
  membership_df <- data.frame("membership" = paste0("Metacell_",na.exclude(seurat.mc.multi@misc$membership)),
                              row.names = colnames(seurat)[!is.na(seurat.mc.multi@misc$membership)])
  
  seurat.mc.multi[[paste0("compactness",opt$Mod1assay)]] <- mc_compactness(cell.membership = membership_df, 
                                                                  sc.obj = seurat,
                                                                  sc.reduction = paste0(opt$Mod1assay,"_diffusion"))
  
  seurat.mc.multi[[paste0("compactness",opt$Mod2assay)]] <- mc_compactness(cell.membership = membership_df, 
                                                                  sc.obj = seurat,
                                                                  sc.reduction = paste0(opt$Mod2assay,"_diffusion"))
  seurat.mc.multi[[paste0("separaton",opt$Mod1assay)]] <- mc_separation(cell.membership = membership_df, 
                                                               sc.obj = seurat,
                                                               sc.reduction =paste0(opt$Mod2assay,"_diffusion"))
  seurat.mc.multi[[paste0("separaton",opt$Mod2assay)]]<- mc_separation(cell.membership = membership_df, 
                                                              sc.obj = seurat,
                                                              sc.reduction = paste0(opt$Mod2assay,"_diffusion"))
  
  seurat.mc.multi$compactness.pca <- mc_compactness(cell.membership = membership_df,sc.obj = seurat, sc.reduction = "pca_diffusion",dims = 1:9)
  seurat.mc.multi$compactness.lsi <- mc_compactness(cell.membership = membership_df,sc.obj = seurat, sc.reduction = "lsi_diffusion",dims = 1:9)
  
  seurat.mc.multi$seperation.pca <-  mc_separation( cell.membership = membership_df,sc.obj = seurat, sc.reduction = "pca_diffusion",dims = 1:9)
  seurat.mc.multi$seperation.lsi <-  mc_separation( cell.membership = membership_df,sc.obj = seurat, sc.reduction = "lsi_diffusion",dims = 1:9)
  meta.data <- seurat.mc.multi@meta.data
  meta.data[[p]] <- NULL
  meta.data$guiding_labels <- p
  meta.data$pct_known <- 100 - pct.to.hide[p] 
  meta.data.all <- rbind(meta.data.all,meta.data)
  
}


ggplot(meta.data.all,aes(x=pct_known,y=seurat_annotations_purity)) + geom_boxplot()

