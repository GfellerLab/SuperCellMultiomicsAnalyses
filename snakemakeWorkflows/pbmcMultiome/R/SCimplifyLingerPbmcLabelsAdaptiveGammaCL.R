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

if (is.null(opt$lingerAggregation)) {
    opt$lingerAggregation <- FALSE
}

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

if(opt$aggregateFragmentfile) {
  outputDirMcFragment <- paste0(opt$outdir,"/aggregated_fragment_file/")
  fragmentFiles <- list()
  fragmentFiles[["ATAC"]] <- GetFragmentData(object = Fragments(seurat)[[1]], slot = "path")
} else {
  fragmentFiles <- NULL
}

SCimplify_for_Seurat_label <- function(label,seurat) {
  seurat.ss <- seurat[,seurat[[opt$label_col]][,1] == label]
  if(ncol(seurat.ss) > 10) {
    N.mc <- floor(sqrt(ncol(seurat.ss)))+1
    g <- ncol(seurat.ss)/N.mc
    DefaultAssay(seurat.ss) <- opt$RNAassay
    seurat.ss <- FindVariableFeatures(seurat.ss) %>% ScaleData() %>% RunPCA()
    DefaultAssay(seurat.ss) <- "RNA"
    seurat.ss <- FindVariableFeatures(seurat.ss) %>% ScaleData() %>% RunPCA()
    
    
    
    if (opt$RNAassay == "SCTransform") {
      rnaAssay = "SCT"
      DefaultAssay(seurat.ss) <- "RNA"
      seurat.ss <- SCTransform(seurat.ss, verbose = FALSE) %>% RunPCA() 
    } else {
      rnaAssay = "RNA"
      seurat.ss <- NormalizeData(seurat.ss, verbose = FALSE) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
    }
    
    DefaultAssay(seurat.ss) <- "ATAC"
    seurat.ss <- RunTFIDF(seurat.ss)
    seurat.ss <- FindTopFeatures(seurat.ss, min.cutoff = "q0")
    seurat.ss <- RunSVD(seurat.ss)

    
    
    seurat.res <- SCimplify_for_Seurat(
      seurat.ss, 
      assay = c(opt$RNAassay,opt$ATACassay),
      k.knn = opt$k.wnn,
      reduction = list("pca", "lsi"), 
      dims = list(opt$RNAcomp, opt$ATACcomp), 
      fragmentFiles = fragmentFiles,
      tmpPath =paste0(opt$outdir,"tmp/"),
      outputDirMcFragment = paste0(getwd(),"/",outputDirMcFragment),
      graph.name = "knn",
      kernel = opt$kernel,
      gamma = g
    )
  } else {
    seurat.res <- seurat.ss
    seurat.res$size <- 1
  }
  return(seurat.res)
}


mc.list <- lapply(X =unique(seurat[[opt$label_col]])[,1], FUN = function(X) {SCimplify_for_Seurat_label(label =X,seurat = seurat)}) 

seurat.mc.multi <- merge(mc.list[[1]],mc.list[-1])

membership <- c()
max.membership <- 0
for(i in 1:length(mc.list)){
  cur.SC <- mc.list[[i]]
  
  cur.membership <- cur.SC@misc$membership
  cur.membership <- cur.membership + max.membership
  membership     <- c(membership, cur.membership)
  max.membership <- max(membership)
}

seurat.mc.multi@misc$membership <- membership

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
  
  
} else {
  DefaultAssay(seurat.mc.multi) <- "ATAC"
  seurat.mc.multi <- RunTFIDF(seurat.mc.multi)
  seurat.mc.multi <- RenameCells(seurat.mc.multi,add.cell.id = "SuperCellMulti")
  atac.data <- seurat.mc.multi[["ATAC"]]@data[linger.prepared.peaks,]
  rna.data <- seurat.mc.multi[["RNA"]]@data[linger.prepared.genes,]
}


#linger.prepared.genes <- linger.prepared.genes[linger.prepared.genes %in% rownames(seurat.mc.multi[["RNA"]])] # not needed when starting from same gene matrix with LingerPrepareCL


rownames(atac.data) <- as.vector(sapply(linger.prepared.peaks,FUN =  function(X){sub("-", ":", X, fixed=TRUE)}))

write.table(rna.data,
            paste0(opt$outdir,"/TG_supercell.tsv"),
            sep = ",",
            quote = F)

write.table(atac.data,paste0(opt$outdir,"/RE_supercell.tsv"),sep = ",",
            quote = F)

saveRDS(seurat.mc.multi, paste0(opt$outdir,"/seurat.multiome.mc.rds"))
