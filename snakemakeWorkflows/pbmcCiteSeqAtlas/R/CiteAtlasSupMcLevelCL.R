library(Seurat)
library(dplyr)
library(getopt)
library(SuperCellMultiomics)


# MetacellExpression <- function(object, pb.method = "aggregate", assays = NULL, features = NULL,
#                                return.seurat = TRUE, group.by = "ident", add.ident = NULL,
#                                layer = "counts", verbose = TRUE, ...)
# {
#   SeuratObject::CheckDots(..., fxns = "CreateSeuratObject")
#   if (!is.null(x = add.ident)) {
#     .Deprecated(msg = "'add.ident' is a deprecated argument, please use the 'group.by' argument instead")
#     group.by <- c("ident", add.ident)
#   }
#   if (!(pb.method %in% c("average", "aggregate"))) {
#     stop("'pb.method' must be either 'average' or 'aggregate'")
#   }
#   object.assays <- .FilterObjects(object = object, classes.keep = c("Assay",
#                                                                     "Assay5"))
#   assays <- assays %||% object.assays
#   if (!all(assays %in% object.assays)) {
#     assays <- assays[assays %in% object.assays]
#     if (length(x = assays) == 0) {
#       stop("None of the requested assays are present in the object")
#     }
#     else {
#       warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
#     }
#   }
#   if (length(x = layer) == 1) {
#     layer <- rep_len(x = layer, length.out = length(x = assays))
#   }
#   else if (length(x = layer) != length(x = assays)) {
#     stop("Number of layers provided does not match number of assays")
#   }
#   data <- FetchData(object = object, vars = rev(x = group.by))
#   data <- data[which(rowSums(x = is.na(x = data)) == 0), ,
#                drop = F]
#   if (nrow(x = data) < ncol(x = object)) {
#     message("Removing cells with NA for 1 or more grouping variables")
#     object <- subset(x = object, cells = rownames(x = data))
#   }
#   for (i in 1:ncol(x = data)) {
#     data[, i] <- as.factor(x = data[, i])
#   }
#   num.levels <- sapply(X = 1:ncol(x = data), FUN = function(i) {
#     length(x = levels(x = data[, i]))
#   })
#   if (any(num.levels == 1)) {
#     message(paste0("The following grouping variables have 1 value and will be ignored: ",
#                    paste0(colnames(x = data)[which(num.levels <= 1)],
#                           collapse = ", ")))
#     group.by <- colnames(x = data)[which(num.levels > 1)]
#     data <- data[, which(num.levels > 1), drop = F]
#   }
#   if (ncol(x = data) == 0) {
#     message("All grouping variables have 1 value only. Computing across all cells.")
#     category.matrix <- matrix(data = 1, nrow = ncol(x = object),
#                               dimnames = list(Cells(x = object), "all"))
#     if (pb.method == "average") {
#       category.matrix <- category.matrix/sum(category.matrix)
#     }
#   }
#   else {
#     category.matrix <- Matrix::sparse.model.matrix(object = as.formula(object = paste0("~0+",
#                                                                                        paste0("data[,", 1:length(x = group.by), "]", collapse = ":"))))
#     print(dim(category.matrix))
#     colsums <- Matrix::colSums(x = category.matrix)
#     category.matrix <- category.matrix[, colsums > 0]
#     colsums <- colsums[colsums > 0]
#     if (pb.method == "average") {
#       category.matrix <- Seurat:::Sweep(x = category.matrix, MARGIN = 2,
#                                         STATS = colsums, FUN = "/")
#     }
#     colnames(x = category.matrix) <- sapply(X = colnames(x = category.matrix),
#                                             FUN = function(name) {
#                                               name <- gsub(pattern = "data\\[, [1-9]*\\]",
#                                                            replacement = "", x = name)
#                                               return(paste0(rev(x = unlist(x = strsplit(x = name,
#                                                                                         split = ":"))), collapse = "_"))
#                                             })
#   }
#   data.return <- list()
#   for (i in 1:length(x = assays)) {
#     data.use <- GetAssayData(object = object, assay = assays[i],
#                              layer = layer[i])
#     features.to.avg <- features %||% rownames(x = data.use)
#     if (inherits(x = features, what = "list")) {
#       features.to.avg <- features[i]
#     }
#     if (IsMatrixEmpty(x = data.use)) {
#       warning("The ", layer[i], " layer for the ", assays[i],
#               " assay is empty. Skipping assay.", immediate. = TRUE,
#               call. = FALSE)
#       next
#     }
#     bad.features <- setdiff(x = features.to.avg, y = rownames(x = data.use))
#     if (length(x = bad.features) > 0) {
#       warning("The following ", length(x = bad.features),
#               " features were not found in the ", assays[i],
#               " assay: ", paste(bad.features, collapse = ", "),
#               call. = FALSE, immediate. = TRUE)
#     }
#     features.assay <- intersect(x = features.to.avg, y = rownames(x = data.use))
#     if (length(x = features.assay) > 0) {
#       data.use <- data.use[features.assay, ]
#     }
#     else {
#       warning("None of the features specified were found in the ",
#               assays[i], " assay.", call. = FALSE, immediate. = TRUE)
#       next
#     }
#     data.return[[i]] <- as.sparse(x = (data.use %*% category.matrix))
#     colnames(data.return[[i]]) <- paste0("Metacell_", c(1:ncol(data.return[[i]])))
#     names(x = data.return)[i] <- assays[[i]]
#   }
#   if (return.seurat) {
#     toRet <- CreateSeuratObject(counts = data.return[[1]],
#                                 project = if (pb.method == "average")
#                                   "Average"
#                                 else "Aggregate", assay = names(x = data.return)[1],
#                                 ...)
#     if (length(x = data.return) > 1) {
#       for (i in 2:length(x = data.return)) {
#         toRet[[names(x = data.return)[i]]] <- CreateAssay5Object(counts = data.return[[i]])
#       }
#     }
#     if (DefaultAssay(object = object) %in% names(x = data.return)) {
#       DefaultAssay(object = toRet) <- DefaultAssay(object = object)
#     }
#     if ("ident" %in% group.by) {
#       first.cells <- c()
#       for (i in 1:ncol(x = category.matrix)) {
#         first.cells <- c(first.cells, Position(x = category.matrix[,
#                                                                    i], f = function(x) {
#                                                                      x > 0
#                                                                    }))
#       }
#       Idents(object = toRet) <- Idents(object = object)[first.cells]
#     }
#     return(toRet)
#   }
#   else {
#     return(data.return)
#   }
# }


set.seed(1234)
spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'gamma', 'g', 1, "numeric", 'gamma to used for metacell identification',
  "unsupIntegration", "u", 0, "logical", "unsupervised integration with STACAS (metacells identificaiton stays supervised)",
  "unsupMetacells", "v", 0, "logical", "unsupervised metacells with SuperCell (integration stays supervised)",
  "sampleName",  "s", 1, "character", "sample name",
  "hide", "e", 1, "numeric", "hide a e proportion of the labels", 
  "dataMetacells", "d", 1, "character", "how to get metacell normalized data: sum_norm (default) or logMean",
  "RNAcomp", "p", 1, "character", "range of components to consider for wnn analysis (eg 1:40 for RNA pca)",
  "ADTcomp", "q", 1, "character", "range of components to consider for wnn analysis (eg 1:50 for ADT pca)",
  "RNAnormalization", "r", 1, "character", "RNA normalisation method for integration  (only 'LogNormalize' at the moment)",
  "RNAmetacells", "m", 1, "character", "RNA normalisation method for metacells identification (default 'LogNormalize' or 'SCT' for SCTransform)",
  "SaveH5adFiles", "f", 1, "logical", "Save H5ad files for the different modalities"
  
), byrow=TRUE, ncol=5)

opt <- list()
#setwd('SuperCellMultiomicsAnalyses/')
# opt$inputSeurat <- "output/pbmcCiteSeqAtlas/pbmc_cite_seq_atlas_filtered.rds"
# opt$RNAnormalization <- "SCT"
# opt$gamma <- 20

opt = getopt(spec)




if(is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "LogNormalize"
}


if(is.null(opt$dataMetacells)) {
  opt$dataMetacells <- "sum_norm"
}

if(is.null(opt$RNAmetacells)) {
  opt$RNAmetacells <- "LogNormalize"
}

if(is.null(opt$SaveH5adFiles)) {
  opt$SaveH5adFiles <- FALSE
  
}


if(is.null(opt$unsupIntegration)) {
  opt$unsupIntegration <- FALSE
  
}

if(is.null(opt$unsupMetacells)) {
  opt$unsupMetacells <- FALSE
  
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


pbmc <- readRDS(opt$inputSeurat)
pbmc$orig.ident <- paste0(pbmc$donor,"_",pbmc$time)
pbmc <- pbmc[,pbmc$celltype.l2 != "Doublet"]
#pbmc <- pbmc[,pbmc$time == 0]  # first test on reference samples
pbmc$celltype.l1.5 <- pbmc$celltype.l1
pbmc$celltype.l1.5[pbmc$celltype.l1 == "other"] <- pbmc$celltype.l2[pbmc$celltype.l1 == "other"]

pbmc$celltype.l1.5[pbmc$celltype.l1 == "other T"] <- pbmc$celltype.l2[pbmc$celltype.l1 == "other T"]
pbmc$label_metacell <- pbmc$celltype.l1.5
if (! is.null(opt$hide)) {
  
  hide.cells <- sample(Cells(pbmc),size = opt$hide*ncol(pbmc))
  pbmc@meta.data[hide.cells,"label_metacell"] <- NA
  if(opt$hide == 1) {
    pbmc$label_metacell <- "unknown" ## Could be in SCimplify: if all NA set to unknown
  }
}

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
    
    pbmcMC <- SCimplify_for_Seurat_v5(seurat = pbmc.smp,
                                      assay = c('RNA','ADT'),
                                      reduction = list("pca", "apca"), 
                                      dims = list(opt$RNAcomp, opt$ADTcomp),
                                      graph.name = "knn",
                                      label = "label_metacell",
                                      kernel = T,
                                      gamma = opt$gamma)
    if (opt$unsupMetacells) {
      pbmcMC <- SCimplify_for_Seurat_v5(seurat = pbmc.smp,
                                        assay = c('RNA','ADT'),
                                        reduction = list("pca", "apca"), 
                                        dims = list(opt$RNAcomp, opt$ADTcomp),
                                        graph.name = "knn",
                                        kernel = T,
                                        gamma = opt$gamma)
    }
    
    if (opt$dataMetacell == "logMean") {
      pbmc.smp <- AddMetaData(pbmc.smp,metadata = pbmcMC@misc$membership[Cells(pbmc.smp)],col.name =  paste0("metacell_g", opt$gamma))
      pbmcMC[["RNA"]]$data <- MetacellExpression(pbmc.smp, 
                                                    assays = "RNA", pb.method = "average",
                                                    group.by = paste0("metacell_g", opt$gamma), 
                                                    layer = "data", 
                                                    return.seurat = F)$RNA
      
      pbmcMC[["ADT"]]$data <- MetacellExpression(pbmc.smp, 
                                                    assays = "ADT", pb.method = "average",
                                                    group.by = paste0("metacell_g", opt$gamma), 
                                                    layer = "data", 
                                                    return.seurat = F)$RNA
      
      if (opt$unsupMetacells) {
        pbmcMC <- SCimplify_for_Seurat_v5(seurat = pbmc.smp,
                                          assay = c('RNA','ADT'),
                                          reduction = list("pca", "apca"), 
                                          dims = list(opt$RNAcomp, opt$ADTcomp),
                                          graph.name = "knn",
                                          kernel = T,
                                          gamma = opt$gamma)
      }
    }
    
    
    
    
  } else {
    pbmc.smp <-  SCTransform(pbmc.smp, verbose = FALSE,conserve.memory = TRUE) %>%
      RunPCA(verbose = FALSE)
    
    pbmcMC <- SCimplify_for_Seurat_v5(seurat = pbmc.smp,
                                      assay = c('SCT','ADT'),
                                      reduction = list("pca", "apca"), 
                                      dims = list(opt$RNAcomp, opt$ADTcomp),
                                      graph.name = "knn",
                                      label = "label_metacell",
                                      kernel = T,
                                      gamma = opt$gamma)
    if (opt$unsupMetacells) {
      pbmcMC <- SCimplify_for_Seurat_v5(seurat = pbmc.smp,
                                        assay = c('RNA','ADT'),
                                        reduction = list("pca", "apca"), 
                                        dims = list(opt$RNAcomp, opt$ADTcomp),
                                        graph.name = "knn",
                                        kernel = T,
                                        gamma = opt$gamma)
    }
    if (opt$dataMetacell == "logMean") {
      DefaultAssay(pbmc.smp) <- "RNA"
      pbmc.smp <- NormalizeData(pbmc.smp)
      pbmc.smp <- AddMetaData(pbmc.smp,metadata = pbmcMC@misc$membership[Cells(pbmc.smp)],col.name =  paste0("metacell_g", opt$gamma))
      pbmcMC[["RNA"]]$data <- MetacellExpression(pbmc.smp, 
                                                    assays = "RNA", pb.method = "average",
                                                    group.by = paste0("metacell_g", opt$gamma), 
                                                    layer = "data", 
                                                    return.seurat = F)$RNA
      
      pbmcMC[["ADT"]]$data <- MetacellExpression(pbmc.smp, 
                                                    assays = "ADT", pb.method = "average",
                                                    group.by = paste0("metacell_g", opt$gamma), 
                                                    layer = "data", 
                                                    return.seurat = F)$ADT
    }
  }
  
  
  remove(pbmc.smp)
  if (length(unique(pbmc$orig.ident)) > 1) {
    pbmc <- pbmc[,pbmc$orig.ident != sampleName]
  }
  gc() 
  
  pbmcMC_supL1.5.list[[sampleName]] <- pbmcMC
  
  
  
}


pbmcMC_supL1.5.list <- lapply(pbmcMC_supL1.5.list, function(x){
  x$label_integration <- x$celltype.l1.5
  if (!is.null(opt$hide)) {
    x[["label_integration"]][,1][x[["label_metacell_with_unknown"]][,1] == "unknown"] <- NA
    if (opt$hide == 1) {
      x$label_integration <- NA
    }
  }
  if (opt$unsupIntegration) {
    x$label_integration <- "None"
  }
  return(x)})


#saveRDS(pbmcMC_supL1.5.list,"pbmcMC_supL1.5.list.rds")


#pbmcMC_supL1.5.list <- readRDS("pbmcMC_supL1.5.list.rds") #identify with prior SCT normalization on RNA seq
library(STACAS)

# integration of ADT data

pbmcMC_supL1.5.list <- lapply(X = pbmcMC_supL1.5.list, FUN = function(x) {
  DefaultAssay(x) <- "ADT"
  x[['SCT']] <- NULL
  VariableFeatures(x) <- rownames(x[["ADT"]])
  x <- RenameCells(x,add.cell.id = unique(x$orig.ident))
  x@misc$new.membership <- paste0(unique(x$orig.ident),"_Metacell_", x@misc$membership)
  names(x@misc$new.membership) <- names(x@misc$membership)
  if (opt$dataMetacells == "logMean") {
    x <- ScaleData(x,features = rownames(x)) %>% RunPCA(reduction.name = 'pca')
    
  } else {
    x <- NormalizeData(x, normalization.method = 'CLR', margin = 2) %>% ScaleData(features = rownames(x)) %>% RunPCA(reduction.name = 'pca')
  }
  return(x)})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = pbmcMC_supL1.5.list)

reference <- which(endsWith(names(pbmcMC_supL1.5.list),suffix = "_0"))

print(reference)

pbmc.anchors <- FindAnchors.STACAS(object.list = pbmcMC_supL1.5.list, 
                                   anchor.features = features,
                                   scale.data = T,
                                   cell.labels = "label_integration",
                                   reference = reference,
                                   dims = opt$ADTcomp)

pbmc.combined <- IntegrateData.STACAS(anchorset = pbmc.anchors,new.assay.name = "integratedADT")

pbmc.combined@misc$membership <- unlist(lapply(pbmcMC_supL1.5.list,FUN = function(x) {x@misc$new.membership}))

Assays(pbmc.combined)


## Integration of RNA data


rnaAssay = "RNA"
pbmcMC_supL1.5.list <- lapply(X = pbmcMC_supL1.5.list, FUN = function(x) { 
  if (opt$dataMetacells == "logMean") {
    DefaultAssay(x) <- "RNA"
  } else {
    DefaultAssay(x) <- "RNA";x <- NormalizeData(x, verbose = FALSE)
  }
  return(x)}
)
features <- SelectIntegrationFeatures(object.list = pbmcMC_supL1.5.list)
pbmcMC_supL1.5.list <- lapply(X = pbmcMC_supL1.5.list, FUN = function(x) {
  x <- ScaleData(x)
  x <- RunPCA(x, features = features, verbose = FALSE)
})




# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features



pbmc.anchors <- FindAnchors.STACAS(object.list = pbmcMC_supL1.5.list, 
                                   anchor.features = features,
                                   cell.labels = "label_integration",
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
pbmc.combined <- RunUMAP(pbmc.combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",min.dist = 0.2)
pbmc.combined <- FindClusters(pbmc.combined, resolution = c(c(5:15)/10), graph.name = "wsnn", algorithm = 3)

saveRDS(pbmc.combined, paste0(opt$outdir,"/seuratCombinedWNN.rds"))




















