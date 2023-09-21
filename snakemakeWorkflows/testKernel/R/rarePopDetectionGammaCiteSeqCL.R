library(Seurat)
library(SuperCellMultiomics)
library(SingleCellExperiment)
library(getopt)
library(future.apply)
library(dplyr)
library(ggplot2)

options(future.globals.maxSize = 10000 * 1024^2)

source("R/functions/rarePopAnalyses.R")


spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : processed seurat object with a dimension reduction computed (.rds)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "cellTypeCol", "c", 1, "character", "colname of cell type analysed",
  "cellType", "t", 1, "character", "name of the cell type for which a rare population is created",
  "RNAcomp", "p", 1, "character", "range of components to consider for metacell identification (eg 1:30 for RNA pca)",
  "ADTcomp", "q", 1, "character", "range of components to consider for metacell identification (eg 1:18 for ADT pca)",
  "k.knn", "k", 1, "numeric", "k for the knn used in the metacell identification",
  "propRarePop", "r", 1, "numeric", "proportion of the rare population that is created",
  "RNAnormalization", 1, "a", "character", "normalisation method for RNA (logNormalize or SCTransform)", 
  "nRep", "n", 1, "numeric", "number of replicates (cell subsampling with different seeds)",
  "nVarGenes", "v", 1, "numeric", "number of variable genes",
  "nWorkers", "w", 1, "numeric", "number of workers to use (default parallel::detectCores()-4)"
), byrow=TRUE, ncol=5)

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code
# test
# setwd("/home/leonard/work/SuperCellMultiomicsAnalyses/")
# opt <- list()
# opt$inputSeurat = "bmcite"
# opt$outdir = "output/testKernel/CITEseq/bmcite/CD56_bright_NK/"
# opt$cellTypeCol = "celltype.l2"
# opt$cellType = "Prog_Mk"
# opt$RNAcomp =  "1:30"
# opt$ADTcomp  = "1:18"
# opt$nVarGenes = 2000
# opt$nWorkers = 4
# opt$k.knn= 30

opt$RNAnormalization =  "logNormalize"


if(is.null(opt$nWorkers)){
  opt$nWorkers = parallel::detectCores()-4
} 

if (is.null(opt$k.knn)) {
  opt$k.knn <- 30
}

if (is.null(opt$nVarGenes)) {
  opt$nVarGenes <- 2000
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

if(is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "logNormalize"
  
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

print(opt)


dir.create(opt$outdir,recursive = T,showWarnings = F)

if(endsWith(opt$inputSeurat,suffix = "rds")) {
  seurat <- readRDS(opt$inputSeurat)
}else{
  seurat <- SeuratData::LoadData(ds = opt$inputSeurat)
}


seurat[[opt$cellTypeCol]][,1] <- gsub(" ","_",seurat[[opt$cellTypeCol]][,1])
seurat[[opt$cellTypeCol]][,1] <- gsub("/","_",seurat[[opt$cellTypeCol]][,1])

DefaultAssay(seurat) <- 'RNA'

if(opt$RNAnormalization == "logNormalize") {
  seurat <- NormalizeData(seurat) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
} else {
  seurat <- SCTransform(seurat) %>% RunPCA()
}


DefaultAssay(seurat) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(seurat) <- rownames(seurat[["ADT"]])
seurat <- NormalizeData(seurat, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')


if(!is.null(opt$propRarePop)) {
  trueProp <-  table(seurat[[opt$cellTypeCol]][,1])[opt$cellType]/ncol(seurat)
  print(trueProp)
  if(opt$propRarePop > trueProp) {
    opt$propRarePop <- NULL
    break
  }
  otherCellNumber <- length(which(seurat[[opt$cellTypeCol]][,1] != opt$cellType))
  
  sizeRarePop <- floor(opt$propRarePop*otherCellNumber/(1-opt$propRarePop))
  print(sizeRarePop)
  
}



## The below is same as plan(multisession, workers=4)
# cl <- parallel::makeCluster(8)
# plan(cluster, workers=cl)


if (!is.null(opt$propRarePop)) {
  plan(multisession, workers=opt$nWorkers)
  
  results <- future_lapply(c(1:opt$nRep), FUN=function(x) {
    tableRes <- data.frame()
    seuratSub <- subsampleCellType(seurat = seurat,
                                   cellType = opt$cellType,
                                   cellTypeCol = opt$cellTypeCol,
                                   n = sizeRarePop,
                                   seed = x,
                                   nfeatures = opt$nVarGenes)
    
    kernelRes <- findRarePop(seuratWithRarePop = seuratSub,
                             cellType = opt$cellType,
                             cellTypeCol = opt$cellTypeCol,
                             kernel = T, 
                             k.knn = opt$k.knn, 
                             assay = c('RNA','ADT'),
                             graph.name = "knn",
                             reduction = list("pca", "apca"),
                             dims = list(opt$RNAcomp,opt$ADTcomp))
    
    tableRes <- rbind(tableRes,c(x,T,unlist(kernelRes),opt$cellType))
    
    standardRes <- findRarePop(seuratWithRarePop = seuratSub,
                               cellType = opt$cellType,
                               cellTypeCol = opt$cellTypeCol,
                               kernel = F, 
                               k.knn = opt$k.knn, 
                               assay = c('RNA','ADT'),
                               graph.name = "knn",
                               reduction = list("pca", "apca"),
                               dims = list(opt$RNAcomp,opt$ADTcomp))
    
    tableRes <- rbind(tableRes,c(x,F,unlist(standardRes),opt$cellType))
    colnames(tableRes) <- c('seed' ,"kernel", 'k' ,'gamma' ,'prop', 'purity',"cellType" )
    return(tableRes)
  }, future.chunk.size=1,future.seed =TRUE)
  
  
  finalTable <- results[[1]]
  
  for (tableRes in results[-1]) {
    finalTable <- rbind(finalTable,tableRes)
  }
} else {
  tableRes <- data.frame()
  x <- 1
  kernelRes <- findRarePop(seuratWithRarePop = seurat,
                           cellType = opt$cellType,
                           cellTypeCol = opt$cellTypeCol,
                           kernel = T, 
                           k.knn = opt$k.knn, 
                           assay = c('RNA','ADT'),
                           graph.name = "knn",
                           reduction = list("pca", "apca"),
                           dims = list(opt$RNAcomp,opt$ADTcomp))
  
  tableRes <- rbind(tableRes,c(x,T,unlist(kernelRes),opt$cellType))
  
  standardRes <- findRarePop(seuratWithRarePop = seurat,
                             cellType = opt$cellType,
                             cellTypeCol = opt$cellTypeCol,
                             kernel = F, 
                             k.knn = opt$k.knn, 
                             assay = c('RNA','ADT'),
                             graph.name = "knn",
                             reduction = list("pca", "apca"), 
                             dims = list(opt$RNAcomp,opt$ADTcomp))
  
  tableRes <- rbind(tableRes,c(x,F,unlist(standardRes),opt$cellType))
  colnames(tableRes) <- c('seed' ,"kernel", 'k' ,'gamma' ,'prop', 'purity',"cellType" )
  finalTable <- tableRes
}



finalTable$seed <- as.numeric(finalTable$seed)
finalTable$gamma <- as.numeric(finalTable$gamma)
finalTable$purity <- as.numeric(finalTable$purity)
finalTable$prop <- as.numeric(finalTable$prop)



finalTable$kernel <- factor(finalTable$kernel ,levels = c(FALSE,TRUE))

pdf(paste0(opt$outdir,"/res.pdf"))

ggplot(finalTable,aes( y = prop, x = kernel,fill = kernel)) + geom_boxplot()  +
  scale_shape(solid = FALSE)  + facet_wrap("~cellType")

ggplot(finalTable,aes( y = purity, x = kernel,fill = kernel)) + geom_boxplot()  +
  scale_shape(solid = FALSE)  + facet_wrap("~cellType")

ggplot(finalTable,aes( y = gamma, x = kernel,fill = kernel)) + geom_boxplot()  +
  scale_shape(solid = FALSE)  + facet_wrap("~cellType")

dev.off()

write.csv(finalTable,paste0(opt$outdir,"/detectionGammaRes.csv"))


