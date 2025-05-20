library(SeuratData)
library(Seurat)
library(Signac)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SuperCellMultiomics)
library(getopt)
library(doParallel)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputRNA_h5ad',  'v', 1, "character", "REQUIRED : path scRNA h5ad preprocessed by LINGER (for annotation and filtered cells/genes)",
  'inputATAC_h5ad',  'w', 1, "character", "REQUIRED : path scATAC h5ad preprocessed by LINGER (for annotation and filtered cells/peaks)",
  'inputSeurat',  'i', 1, "character", "REQUIRED : path preprocessed seurat multiomic rds object",
  'label_col', "l", 1, "character", "name of label column to guide metacell identification",
  'lingerAggregation',"n", 0, "logical", "Use linger aggregartion (average RNA lognorm, average ATAC log1p_raw).",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "RNAcomp", "p", 1, "character", "range of RNA components to consider for metacell identification  (eg 1:50 for RNA pca)",
  "ATACcomp", "q", 1, "character", "range of ATAC components to consider for metacell identification (eg 2:50 for ATAC pca)",
  "RNAassay", "r", 1, "character", "RNA assay name with computed pca (default RNA)",
  "ATACassay", "a", 1, "character", "ATAC assay name with computed lsi (default ATAC)" , 
  "k.wnn", "k", 1, "numeric", "k for the knn used in the wnn used for metacell identification",
  "gamma", "g", 1, "numeric", "gamma used for metacell identificaiton",
  "kernel", "e", 1, "logical", "whether to use a kernel",
  "randomMetacells", "d", 1, "logical","construct random metacell at the specified gamma",
  "memberships", "c", 1 , "character", "memberships already computed (e.g. with SEACells) csv table: 1st column cell name, 2nd column metacell name",
  "prefixMembership", "x", 1, "character", "metacell name prefix memebership to disccard (e.g. SEACell-)",
  "inputSeuratMetacell", "m", 1, "character", "seurat metacell object if available to rescale directly",
  "aggregateFragmentfile", "f", 1, "character", 'whether to aggregate fragment file or not'
  
), byrow=TRUE, ncol=5)

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code
# test
# opt <- list()
# opt$aggregateFragmentfile <-T
# opt$inputSeurat <- "output/correlationAnalyzis/pbmcMultiome/singlecells_analysis/seuratWNN.rds"
# opt$outdir <- "output/correlationAnalyzis/pbmcMultiome/randomMetacells/"
# opt$RNAcomp <- "1:50"
# opt$ATACcomp <- "2:50"
# opt$kernel <- TRUE
# opt$aggregateFragmentfile <- TRUE
# opt$gamma <- 50
# opt$k.wnn <- 30
# opt$randomMetacells = T
#outputDirMcFragment <- "~/work/SuperCellMultiomicsAnalyses/input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
#fragmentFiles <- list()
#fragmentFiles[["ATAC"]] <- "~/work/SuperCellMultiomicsAnalyses/input/pbmcMultiome/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

if (is.null(opt$RNAassay)) {
  opt$RNAassay <- "RNA"
}

if (is.null(opt$ATACassay)) {
  opt$ATACassay <- "ATAC"
}

if(is.null(opt$k.wnn)) {
  opt$k.wnn <- 30
}

if (is.null(opt$lingerAggregation)) {
  opt$lingerAggregation <- FALSE
}


if(is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "logNormalize"
  
}

if (!is.null(opt$RNAcomp)) {
  ci <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][2])
  opt$RNAcomp <- c(ci:cf)
}

if (!is.null(opt$ATACcomp)) {
  ci <- as.numeric(strsplit(opt$ATACcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$ATACcomp,split = ":")[[1]][2])
  opt$ATACcomp <- c(ci:cf)
}

if (is.null(opt$outdir)) {
  opt$outdir <- "./"
}

if (is.null(opt$randomMetacells)) {
  opt$randomMetacells <- F
}

print(opt)

dir.create(opt$outdir,recursive = T,showWarnings = F)


linger.prepared.obs <- anndata::read_h5ad(opt$inputRNA_h5ad,backed = "r")$obs
linger.prepared.genes <- anndata::read_h5ad(opt$inputRNA_h5ad,backed = "r")$var_names
linger.prepared.peaks <- anndata::read_h5ad(opt$inputATAC_h5ad,backed = "r")$var$gene_ids
linger.prepared.peaks <- as.vector(sapply(linger.prepared.peaks,FUN =  function(X){sub(":", "-", X, fixed=TRUE)}))



seurat <- readRDS(opt$inputSeurat)
seurat <- seurat[,linger.prepared.obs$barcode] # keep only cell used by LINGER
rownames(linger.prepared.obs) <- linger.prepared.obs$barcode
seurat$label <- linger.prepared.obs[colnames(seurat),"label"]

## Seurat preprocessing after linger features filtering (not needed)

# if (opt$RNAnormalization == "SCTransform") {
#   rnaAssay = "SCT"
#   DefaultAssay(seurat) <- "RNA"
#   seurat <- SCTransform(seurat, verbose = FALSE) %>% RunPCA(reduction.name = "pca.linger.genes") 
# } else {
#   rnaAssay = "RNA"
#   seurat <- NormalizeData(seurat, verbose = FALSE) %>% FindVariableFeatures(pbmc,nFeature = opt$nVarGenes) %>% ScaleData(pbmc) %>% RunPCA(reduction.name = "pca.linger.genes")
# }
# 
# DefaultAssay(seurat) <- "ATAC"
# seurat <- RunTFIDF(seurat)
# seurat <- FindTopFeatures(seurat, min.cutoff = "q0")
# seurat <- RunSVD(seurat,reduction.name = "lsi.linger.peaks")


  
if (is.null(opt$inputSeuratMetacell)) {
  if(opt$aggregateFragmentfile) {
    outputDirMcFragment <- paste0(opt$outdir,"/aggregated_fragment_file/")
    fragmentFiles <- list()
    fragmentFiles[["ATAC"]] <- GetFragmentData(object = Fragments(seurat)[[1]], slot = "path")
  } else {
    fragmentFiles <- NULL
  }
  
  print(fragmentFiles)
  
  if (opt$randomMetacells) {  
    print("constructing random metacell")
    randomMC <- 1:floor(ncol(seurat)/opt$gamma)
    randomMemberships <-sample(c(randomMC, sample(randomMC, ncol(seurat)-length(randomMC), replace=TRUE)))
    names(randomMemberships) <- colnames(seurat)
    seurat.mc.multi <- SCimplify_for_Seurat(seurat = seurat,
                                            membership=randomMemberships,
                                            fragmentFiles = fragmentFiles,
                                            tmpPath =paste0(opt$outdir,"tmp/"),
                                            outputDirMcFragment = paste0(getwd(),"/",outputDirMcFragment))
  }
  if (!is.null(opt$memberships)) {  
    print("aggregating data with the given memberships")
    membershipsTable <- read.csv(opt$memberships)
    if (is(membershipsTable[,2])[1] == "numeric") { #for MetaCell memberships with discarded outliers with a negative values
      membershipsTable[membershipsTable$membership < 0,2] <- NA
      memberships <- membershipsTable[,2]
    } else { # for SEACells memberships with the names of the archetypes in characters
      memberships <- as.numeric(plyr::mapvalues(x = membershipsTable[,2],from = unique(membershipsTable[,2]), to = c(1:length(unique(membershipsTable[,2])))))  
    }
    
    # 
    
    names(memberships) <- membershipsTable[,1]
    
    #AggregateExpression in SCimplify_for_Seurat remove cells with NA (MetaCell outliers) before aggregating expression
    seurat.mc.multi <- SCimplify_for_Seurat(seurat = seurat,
                                            membership=memberships,
                                            fragmentFiles = fragmentFiles,
                                            tmpPath =paste0(opt$outdir,"tmp/"),
                                            outputDirMcFragment = paste0(getwd(),"/",outputDirMcFragment))
  }
  
  if (!is.null(opt$RNAcomp)|!is.null(opt$ATACcomp)) {
    
    if (!is.null(opt$RNAcomp)&!is.null(opt$ATACcomp)) {
      
      seurat.mc.multi <- SCimplify_for_Seurat(
        seurat, 
        label = opt$label_col,
        assay = c(opt$RNAassay,opt$ATACassay),
        k.knn = opt$k.wnn,
        reduction = list("pca", "lsi"), 
        dims = list(opt$RNAcomp, opt$ATACcomp), 
        fragmentFiles = fragmentFiles,
        tmpPath =paste0(opt$outdir,"tmp/"),
        outputDirMcFragment = paste0(getwd(),"/",outputDirMcFragment),
        graph.name = "knn",
        kernel = opt$kernel,
        gamma = opt$gamma
      )
      } else {
        if (is.null(opt$RNAcomp)&!is.null(opt$ATACcomp)) {
          print("identifying metacells on ATAC modality")
          seurat.mc.multi <- SCimplify_for_Seurat(
            seurat, 
            label = opt$label_col,
            assay = c(opt$ATACassay),
            k.knn = opt$k.wnn,
            reduction = list("lsi"), 
            dims = list(opt$ATACcomp), 
            fragmentFiles = fragmentFiles,
            tmpPath =paste0(opt$outdir,"tmp/"),
            outputDirMcFragment = paste0(getwd(),"/",outputDirMcFragment),
            kernel = opt$kernel,
            gamma = opt$gamma
          )
        }
        
        if (!is.null(opt$RNAcomp)&is.null(opt$ATACcomp)) {
          print("identifying metacells on RNA modality")
          seurat.mc.multi <- SCimplify_for_Seurat(
            seurat, 
            assay = c(opt$RNAassay),
            k.knn = opt$k.wnn,
            reduction = list("pca"), 
            dims = list(opt$RNAcomp), 
            fragmentFiles = fragmentFiles,
            label = opt$label_col,
            tmpPath =paste0(opt$outdir,"tmp/"),
            outputDirMcFragment = paste0(getwd(),"/",outputDirMcFragment),
            kernel = opt$kernel,
            gamma = opt$gamma)
        }
      }
  } 
} else {
  seurat.mc <- readRDS(opt$inputSeuratMetacell)
  
  
  seurat.mc.multi <- SCimplify_for_Seurat(seurat = seurat, seurat.mc = seurat.mc, gamma = opt$gamma)
  
}


DefaultAssay(seurat.mc.multi) <- "RNA"
seurat.mc.multi <- NormalizeData(seurat.mc.multi)

if (opt$lingerAggregation) {
  print("use linger aggregation strategy")
  DefaultAssay(seurat) <- "RNA"
  seurat <- NormalizeData(seurat)
  #seurat$metacells <- seurat.mc.multi@misc$membership[colnames(seurat)]
  rna.data <- supercell_GE(seurat[["RNA"]]@data,groups = seurat.mc.multi@misc$membership[colnames(seurat)],mode = "average")
  rna.data <- rna.data[linger.prepared.genes,]
  colnames(rna.data) <- paste0("SuperCellMulti_",c(1:ncol(rna.data)))
  atac.log1p.sc <- log1p(seurat[["ATAC"]]@counts)[linger.prepared.peaks,]
  atac.data <- supercell_GE(atac.log1p.sc,groups = seurat.mc.multi@misc$membership[colnames(seurat)],mode = "average")
  atac.data[atac.data>100] <- 100
  colnames(atac.data) <- paste0("SuperCellMulti_",c(1:ncol(atac.data)))
  seurat.mc.multi <- RenameCells(seurat.mc.multi,add.cell.id = "SuperCellMulti")
  
  
} else { # TF-IDF certainly inappropriate
  DefaultAssay(seurat.mc.multi) <- "ATAC"
  seurat.mc.multi <- RunTFIDF(seurat.mc.multi)
  seurat.mc.multi <- RenameCells(seurat.mc.multi,add.cell.id = "SuperCellMulti")
  atac.data <- seurat.mc.multi[["ATAC"]]@data[linger.prepared.peaks,]
  rna.data <- seurat.mc.multi[["RNA"]]@data[linger.prepared.genes,]
}




# seurat.mc.multi <- RunTFIDF(seurat.mc.multi)
seurat.mc.multi <- RenameCells(seurat.mc.multi,add.cell.id = "SuperCellMulti")
#linger.prepared.genes <- linger.prepared.genes[linger.prepared.genes %in% rownames(seurat.mc.multi[["RNA"]])] # not needed when starting from same gene matrix with LingerPrepareCL


rownames(atac.data) <- as.vector(sapply(linger.prepared.peaks,FUN =  function(X){sub("-", ":", X, fixed=TRUE)}))

write.table(seurat.mc.multi[["RNA"]]@data[linger.prepared.genes,],
            paste0(opt$outdir,"/TG_supercell.tsv"),
            sep = ",",
            quote = F)

write.table(atac.data,paste0(opt$outdir,"/RE_supercell.tsv"),sep = ",",
            quote = F)

saveRDS(seurat.mc.multi, paste0(opt$outdir,"/seurat.multiome.mc.rds"))
