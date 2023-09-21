library(Seurat)
library(SuperCellMultiomics)
library(SingleCellExperiment)
library(getopt)
library(future.apply)


source("R/functions/rarePopAnalyses.R")


spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : processed seurat object with a dimension reduction computed (.rds)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "cellTypeCol", "c", 1, "character", "colname of cell type analysed",
  "cellType", "t", 1, "character", "name of the cell type for which a rare population is created",
  "components", "p", 1, "character", "range of components to consider for metacell identification (eg 1:30 for RNA pca or 2:30 for ATAC lsi)",
  "k.knn", "k", 1, "numeric", "k for the knn used in the metacell identification",
  "propRarePop", "r", 1, "numeric", "proportion of the rare population that is created",
  "nRep", "n", 1, "numeric", "number of replicates (cell subsampling with different seeds)",
  "nVarGenes", "v", 1, "numeric", "number of variable genes",
  "gamma", "g", 1, "numeric", "gamma for metacell identification",
  "kernel", "d", 1, "logical", "wether to use a kernel or not for metacell identification",
  "nWorkers", "w", 1, "numeric", "number of workers to use (default parallel::detectCores()-4)"
), byrow=TRUE, ncol=5)

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code
# test
# setwd("/home/leonard/work/SuperCellMultiomicsAnalyses/")
# opt <- list()
# opt$outdir <- "./output/testKernel/sc_mixology/H1975"
# opt$cellType <- "H1975"
# opt$inputSeurat <- "input/sce_sc_10x_5cl_qc.rds"
# opt$propRarePop <- 0.005
# opt$cellTypeCol <- "cell_line"
# opt$k.knn = 30
# opt$nRep = 2
# opt$components = "1:30"
# opt$nVarGenes = 2000
# opt$gamma <- 20
# opt$nWorkers = 5
# opt$kernel = T

if (is.null(opt$k.knn)) {
  opt$k.knn <- 30
}

if (is.null(opt$gamma)) {
  opt$gamma <- 20
}

if (is.null(opt$kernel)) {
  opt$kernel <- T
}

if (is.null(opt$nVarGenes)) {
  opt$nVarGenes <- 2000
}

if (is.null(opt$components)) {
  opt$components <- c(1:30)
} else {
  ci <- as.numeric(strsplit(opt$components,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$components,split = ":")[[1]][2])
  opt$components <- c(ci:cf)
}

if(is.null(opt$nWorkers)){
  opt$nWorkers = parallel::detectCores()-4
} 


args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$inputSeurat)) {
  cat("Experiment for detection of artificial rare cell population")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}







if (is.null(opt$outdir)) {
  opt$outdir = "./"
}

dir.create(opt$outdir,recursive = T,showWarnings = F)

seurat <- readRDS(opt$inputSeurat)

seurat <- FindVariableFeatures(seurat,nfeatures = opt$nVarGenes)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)


res0 <- data.frame()




plan(multisession, workers=opt$nWorkers)


results <- future_lapply(c(1:opt$nRep), FUN=function(x) {
  notDetected <- T
  sizePop <- 0
  res <- data.frame(seed = x,
                    kernel = opt$kernel,
                    rarePopSize = NA,
                    rarePopProp = NA)
  
  while(notDetected) {
    sizePop <- sizePop + 1
    print(paste0("Subsampling to ", sizePop, " cells"))
    seuratWithRarePop  <- subsampleCellType(seurat = seurat,
                                            cellType = opt$cellType,
                                            cellTypeCol = opt$cellTypeCol,
                                            n = sizePop,
                                            seed = x,
                                            nfeatures = opt$nVarGenes)
    
    
    rarePopSC  <-  SCimplify_for_Seurat(seurat = seuratWithRarePop,
                                        gamma = opt$gamma,
                                        kernel = opt$kernel, 
                                        k.knn = opt$k.knn, 
                                        dims = list(opt$components),
                                        return.seurat = F)
    
    rarePopSC[[opt$cellTypeCol]] <- supercell_assign(clusters = seuratWithRarePop[[opt$cellTypeCol]][,1],
                                                     supercell_membership = rarePopSC$membership,
                                                     method = "absolute")
    
    rarePopSC[["purity"]] <- supercell_purity(clusters = seuratWithRarePop[[opt$cellTypeCol]][,1], 
                                              supercell_membership  = rarePopSC$membership)

    p <- rarePopSC$purity[rarePopSC[[opt$cellTypeCol]] == opt$cellType]
    print(length(which(p > 0.8))) 
    if (length(which(p > 0.8))>0) {
      res$rarePopSize <- sizePop
      res$rarePopProp<- sizePop/ncol(seuratWithRarePop)
      notDetected <- FALSE
    }
  }  
  return(res)
}, future.chunk.size=1,future.seed =TRUE)


finalTable <- results[[1]]

for (tableRes in results[-1]) {
  finalTable <- rbind(finalTable,tableRes)
}



write.csv(finalTable,paste0(opt$outdir,"/detectionSizeRes.csv"))







