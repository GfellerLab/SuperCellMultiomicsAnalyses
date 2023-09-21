library(SeuratData)
library(Seurat)
library(Signac)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SuperCellMultiomics)
library(getopt)
library(doParallel)
library(chromVAR) #
library(JASPAR2020) #
library(TFBSTools) #
library(motifmatchr) #
library(ggplot2)

source("R/functions/innerNormalizedVariance.R")
source("R/functions/metacellCompactnessSeparation.R")


spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds) or dataset name (from SeuratData package)",
  'geneList', "g", 1, "character", "REQUIRED: gene list to test correlation between gene body accessibility and gene expression (.txt)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'useSize', 's', 1, "logical", "whether to use (metacell) size to compute correlation",
  'chromVAR', "c", 1, "logical", "whether to do a chromVar analyzis or not",
  'nWorkers', "w", 1, "numeric", "number of core to use for chromVAR analysis",
  "compactnessSeparation","p", 1, "logical", "whether to compute compactness and separation using palantir/seacells or not",
  "pythonSeacellEnv", "e", 1, "character", "python path for seacell env to compute compactness/separation",
  "singleCellSeurat", "l", 1, "character", "single cell seurat object to compute compactness/separation",
  "RNAcomp", "n", 1, "character", "range of RNA components to consider for metacell identification  (eg 1:50 for RNA pca)",
  "ATACcomp", "q", 1, "character", "range of ATAC components to consider for metacell identification (eg 2:50 for ATAC pca)"
  
), byrow=TRUE, ncol=5)



#opt <- list()
# opt$inputSeurat <- "output/correlationAnalyzis/pbmcMultiome/metacells/g50/seurat.multiome.mc.rds"
# opt$geneList <- "output/correlationAnalyzis/pbmcMultiome/selectedGenes.txt"
opt = getopt(spec)

if(is.null(opt$nWorkers)){
  opt$nWorkers = parallel::detectCores()-4
} 

if(is.null(opt$compactnessSeparation)){
  opt$compactnessSeparation = F
} 



print(opt)

seurat <- readRDS(opt$inputSeurat)
genes <- read.table(opt$geneList)$x

DefaultAssay(seurat) <- "RNA"
seurat  <- NormalizeData(seurat) 

DefaultAssay(seurat) <- "ATAC"
gene.activities <- GeneActivity(seurat,features =genes )


# add gene activities as a new assay
seurat[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(seurat) <- "ACTIVITY"
seurat <- NormalizeData(seurat)
#seurat <- ScaleData(seurat, features = rownames(seurat))
corrTable <- supercell_FeatureFeaturePlot_Seurat(seurat.mc = seurat,
                                                 feature_x = rownames(seurat[["ACTIVITY"]]),
                                                 feature_y = rownames(seurat[["ACTIVITY"]]),
                                                 assays = c("ACTIVITY",'RNA'),
                                                 method = "spearman",
                                                 is.normalized= T,
                                                 use.size = opt$useSize,
                                                 cluster = "orig.ident", plot = F)

if (opt$chromVAR) {
  library(BiocParallel)
  register(MulticoreParam(opt$nWorkers))
  
  DefaultAssay(seurat) <- "ATAC"
  
  
  # Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
  pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
  motif.matrix <- CreateMotifMatrix(features = granges(seurat), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
  motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
  seurat <- SetAssayData(seurat, assay = 'ATAC', slot = 'motifs', new.data = motif.object)
  
  # Note that this step can take 30-60 minutes 
  seurat <- RunChromVAR(
    object = seurat,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  
}

if (opt$compactnessSeparation) {
  print("computing compactness and separation & inner normalized variance")
  
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
  
  seurat.sc <- readRDS(opt$singleCellSeurat)
  library(reticulate)
  use_python(opt$pythonSeacellEnv)
  source_python("config/SuperCellMultiomics/inst/seacellsBenchmarkMetrics.py")
  
  
  seurat <- metacellCompactnessSeparation(sc.seurat = seurat.sc,SC_seurat = seurat,dims = opt$RNAcomp)
  seurat <- metacellCompactnessSeparation(sc.seurat = seurat.sc,SC_seurat = seurat,preprocessing_method = "lsi", opt$ATACcomp)
  ## Inner normalized variance
  
  seurat$innerNormVar <- computeInnerNormVar(seurat = seurat.sc, memberships = seurat@misc$membership)
  
  
}



## Cell cycle analysis
DefaultAssay(seurat) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

if (!opt$useSize) {
  seurat$size <- 1
}
pdf(paste0(opt$outdir,"/cellCyclePLot.pdf"))
plot(ggplot2::ggplot(seurat@meta.data,aes(x=size,y=G2M.Score,color = Phase)) + geom_point())
dev.off()


## TF expression / TF motif accessibility correlation analyzis
DefaultAssay(seurat) <- "ATAC"
detectedMotifs <- rownames(seurat@assays$chromvar) #motif
relatedTFs <- ConvertMotifID(seurat, id = detectedMotifs) #tf

feature_x <- names(rowSums(seurat@assays$RNA[rownames(seurat@assays$RNA) %in% relatedTFs,]) > 0)[rowSums(seurat@assays$RNA[rownames(seurat@assays$RNA) %in% relatedTFs,]) > 0] # expressed tfs

feature_y <- ConvertMotifID(seurat, name = feature_x) # related motifs

DefaultAssay(seurat) <- "RNA"


# DefaultAssay(seurat.mc.multi) <- "ACTIVITY"
# seurat.mc.multi <- NormalizeData(seurat.mc.multi)

chromVarCorrTable <- supercell_FeatureFeaturePlot_Seurat(seurat.mc = seurat,
                                                             is.normalized = T,
                                                             cluster = "seurat_annotations",
                                                             assays = c("RNA","chromvar"),
                                                             method = "spearman",
                                                             feature_x = feature_x,
                                                             feature_y = feature_y,
                                                             use.size = opt$useSize,
                                                             plot = F)





saveRDS(seurat,paste0(opt$outdir,"/seurat.multiome.mc.activities.rds"))


write.csv(corrTable,paste0(opt$outdir,"/corrTable.csv"))

write.csv(chromVarCorrTable,paste0(opt$outdir,"/chromVarCorrTable.csv"))