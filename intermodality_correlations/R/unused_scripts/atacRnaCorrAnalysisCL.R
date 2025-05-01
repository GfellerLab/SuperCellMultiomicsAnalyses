library(Seurat)
library(Signac)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SuperCellMultiomics)
library(getopt)
library(doParallel)
library(chromVAR) #
# library(JASPAR2020) #
library(TFBSTools) #
library(motifmatchr) #
library(ggplot2)
library(MetacellAnalysisToolkit)
# source("R/functions/innerNormalizedVariance.R")
# source("R/functions/metacellCompactnessSeparation.R")


spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds) or dataset name (from SeuratData package)",
  'geneList', "g", 1, "character", "REQUIRED: gene list to test correlation between gene body accessibility and gene expression (.txt)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'useSize', 's', 1, "logical", "whether to use (metacell) size to compute correlation",
  'pwmChromVAR', "c", 1, "character", "pwm motif matrix for chromvar (if not no chromvar analysis",
  'macsPeaks', 'm', 0, "logical", "To be used with -c, add -m to call macs peak on each single-cell seurat cluster (res 0.8), these peaks will be use for chromvar analysis",
  'nWorkers', "w", 1, "numeric", "number of core to use for chromVAR analysis",
  "compactnessSeparation","p", 0, "logical", "add -p to compute compactness and separation using palantir/seacells or not",
  "pythonSeacellEnv", "e", 1, "character", "python path for seacell env to compute compactness/separation",
  "singleCellSeurat", "l", 1, "character", "single cell seurat object to compute compactness/separation",
  "RNAcomp", "n", 1, "character", "range of RNA components to consider for metacell identification  (eg 1:50 for RNA pca)",
  "ATACcomp", "q", 1, "character", "range of ATAC components to consider for metacell identification (eg 2:50 for ATAC pca)"
), byrow=TRUE, ncol=5)


opt = getopt(spec)

# opt <- list()
# #setwd("~/work/SuperCellMultiomicsAnalyses/")
# opt$inputSeurat <- "output/pbmcMultiome/SuperCellMulti/g10/seurat.multiome.mc.rds"
# opt$singleCellSeurat <- "output/pbmcMultiome/singlecells_analysis/seuratWNN_with_macs2_peaks.rds"
# opt$geneList <- "output/correlationAnalyzis/pbmcMultiome/selectedGenes.txt"
# opt$pwmChromVAR <- "input/JASPAR_2024_human_motifs/pwm.rds"
# opt$nWorkers <- 4
# opt$compactnessSeparation <- T
# opt$macsPeaks <- T
# opt$pythonSeacellEnv <- "/work/FAC/FBM/LLB/dgfeller/scrnaseq/lherault/bin/miniconda3/envs/seaMetaCells"
# opt$RNAcomp <- "1:40"
# opt$ATACComp <- "2:40"
# opt$useSize <- T



if(is.null(opt$macsPeaks)){
  opt$macsPeaks <- F
} 

if(is.null(opt$nWorkers)){
  opt$nWorkers = parallel::detectCores()-4
} 

if(is.null(opt$compactnessSeparation)){
  opt$compactnessSeparation = F
} 



print(opt)

seurat.mc <- readRDS(opt$inputSeurat)
# 
DefaultAssay(seurat.mc) <- "RNA"
seurat.mc  <- NormalizeData(seurat.mc)
# 
DefaultAssay(seurat.mc) <- "ATAC"

if(!"ACTIVITY" %in% Assays(seurat.mc)) {
  genes <- read.table(opt$geneList)$x
  gene.activities <- GeneActivity(seurat.mc,features =genes )
  seurat.mc[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
}


# add gene activities as a new assay

# normalize gene activities
DefaultAssay(seurat.mc) <- "ACTIVITY"
seurat.mc <- NormalizeData(seurat.mc)
#seurat.mc <- ScaleData(seurat.mc, features = rownames(seurat.mc))
corrTable <- supercell_FeatureFeaturePlot_Seurat(seurat.mc = seurat.mc,
                                                 feature_x = rownames(seurat.mc[["ACTIVITY"]]),
                                                 feature_y = rownames(seurat.mc[["ACTIVITY"]]),
                                                 assays = c("ACTIVITY",'RNA'),
                                                 method = "spearman",
                                                 is.normalized= T,
                                                 use.size = opt$useSize,
                                                 cluster = "orig.ident", plot = F)

corrTablePearson <- supercell_FeatureFeaturePlot_Seurat(seurat.mc = seurat.mc,
                                                        feature_x = rownames(seurat.mc[["ACTIVITY"]]),
                                                        feature_y = rownames(seurat.mc[["ACTIVITY"]]),
                                                        assays = c("ACTIVITY",'RNA'),
                                                        method = "pearson",
                                                        is.normalized= T,
                                                        use.size = opt$useSize,
                                                        cluster = "orig.ident", plot = F)

if (!is.null(opt$singleCellSeurat)) { # directly single cells in input
  seurat.sc <- readRDS(opt$singleCellSeurat) # not the most efficient with two sc objects load but avoid writing new code
}

if (!is.null(opt$pwmChromVAR)) {
  
  if (opt$macsPeaks) { 
    peakAssay <- "peaks"
    if (!is.null(opt$singleCellSeurat)) {
      
      
      seurat.sc$metacell <- seurat.mc@misc$membership
      seurat.mc.macs2 <- Seurat::AggregateExpression(seurat.sc, 
                                                     assays = 'peaks', 
                                                     group.by = "metacell", 
                                                     slot = "counts", 
                                                     return.seurat = T)
      
      seurat.mc[['peaks']] <- Signac::as.ChromatinAssay(x = seurat.mc.macs2[['peaks']], 
                                                        genome = unique(seurat.sc@assays[["peaks"]]@seqinfo@genome)[1], 
                                                        ranges = StringToGRanges(rownames(seurat.mc.macs2[["peaks"]]), 
                                                                                 sep =  c("-", "-")), 
                                                        fragments = Fragments(seurat.sc[["peaks"]]), 
                                                        annotation = Annotation(seurat.sc[["peaks"]]))
      DefaultAssay(seurat.mc) <- peakAssay
      remove(seurat.mc.macs2)
      gc()
      
      
      ## add new aggregated peaks to metacells
      
      
    } else {  
      DefaultAssay(seurat.mc) <- peakAssay
    } 
  } else {
    peakAssay <- "ATAC"
    DefaultAssay(seurat.mc) <- peakAssay
  }
  
  
  library(BiocParallel)
  register(MulticoreParam(opt$nWorkers))
  
  
  
  # Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
  #pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
  print('reading pwm')
  pwm_set <- readRDS(opt$pwmChromVAR )
  print("creating motif matrix")
  motif.matrix <- CreateMotifMatrix(features = granges(seurat.mc), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
  print("creating motif object")
  motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
  print("setting assay data")
  seurat.mc <- SetAssayData(seurat.mc, assay = peakAssay, slot = 'motifs', new.data = motif.object)
  
  # Note that this step can take 30-60 minutes 
  seurat.mc <- RunChromVAR(
    object = seurat.mc,
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
  
  library(reticulate)
  use_python(opt$pythonSeacellEnv)
  #source_python("config/SuperCellMultiomics/inst/seacellsBenchmarkMetrics.py")
  membership.df <- data.frame('membership' = seurat.mc@misc$membership)
  # library(MetacellAnalysisToolkit)
  # excluded.genes <- rownames(pbmc[['RNA']])[startsWith(rownames(pbmc[['RNA']]),prefix = "MT-")] #only work for human
  # 
  # excluded.genes <- c(excluded.genes,rownames(pbmc[['RNA']])[startsWith(rownames(pbmc[['RNA']]),prefix = "RPL")]) #only work for human
  # 
  # excluded.genes <- c(excluded.genes,rownames(pbmc[['RNA']])[startsWith(rownames(pbmc[['RNA']]),prefix = "RPS")]) #only work for human
  # 
  # excluded.genes <- c(excluded.genes,c("ACSM3", "ANP32B", "APOE", "AURKA", "B2M", "BIRC5", "BTG2", "CALM1", "CD63", "CD69", "CDK4","CENPF", "CENPU", "CENPW", "CH17-373J23.1", "CKS1B", "CKS2", "COX4I1", "CXCR4", "DNAJB1",
  #                                      "DONSON", "DUSP1", "DUT", "EEF1A1", "EEF1B2", "EIF3E", "EMP3", "FKBP4", "FOS", "FOSB", "FTH1",
  #                                      "G0S2", "GGH", "GLTSCR2", "GMNN", "GNB2L1", "GPR183", "H2AFZ", "H3F3B", "HBM", "HIST1H1C",
  #                                      "HIST1H2AC", "HIST1H2BG", "HIST1H4C", "HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB",
  #                                      "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-E", "HLA-F", "HMGA1",
  #                                      "HMGB1", "HMGB2", "HMGB3", "HMGN2", "HNRNPAB", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B",
  #                                      "HSPA6", "HSPD1", "HSPE1", "HSPH1", "ID2", "IER2", "IGHA1", "IGHA2", "IGHD", "IGHG1", "IGHG2",
  #                                      "IGHG3", "IGHG4", "IGHM", "IGKC", "IGKV1-12", "IGKV1-39", "IGKV1-5", "IGKV3-15", "IGKV4-1",
  #                                      "IGLC2", "IGLC3", "IGLC6", "IGLC7", "IGLL1", "IGLL5", "IGLV2-34", "JUN", "JUNB", "KIAA0101",
  #                                      "LEPROTL1", "LGALS1", "LINC01206", "LTB", "MCM3", "MCM4", "MCM7", "MKI67", "MT2A", "MYL12A",
  #                                      "MYL6", "NASP", "NFKBIA", "NUSAP1", "PA2G4", "PCNA", "PDLIM1", "PLK3", "PPP1R15A", "PTMA",
  #                                      "PTTG1", "RAN", "RANBP1", "RGCC", "RGS1", "RGS2", "RGS3", "RP11-1143G9.4", "RP11-160E2.6",
  #                                      "RP11-53B5.1", "RP11-620J15.3", "RP5-1025A1.3", "RP5-1171I10.5", "RPS10", "RPS10-NUDT3", "RPS11",
  #                                      "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17", "RPS18", "RPS19", "RPS19BP1",
  #                                      "RPS2", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS27L",
  #                                      "RPS28", "RPS29", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS4Y2", "RPS5", "RPS6", "RPS6KA1",
  #                                      "RPS6KA2", "RPS6KA2-AS1", "RPS6KA3", "RPS6KA4", "RPS6KA5", "RPS6KA6", "RPS6KB1", "RPS6KB2",
  #                                      "RPS6KC1", "RPS6KL1", "RPS7", "RPS8", "RPS9", "RPSA", "RRM2", "SMC4", "SRGN", "SRSF7", "STMN1",
  #                                      "TK1", "TMSB4X", "TOP2A", "TPX2", "TSC22D3", "TUBA1A", "TUBA1B", "TUBB", "TUBB4B", "TXN", "TYMS",
  #                                      "UBA52", "UBC", "UBE2C", "UHRF1", "YBX1", "YPEL5", "ZFP36", "ZWINT")) #LATERAL genes MC2
  # 
  # excluded.genes <- c(excluded.genes,c("CCL3", "CCL4", "CCL5", "CXCL8", "DUSP1", "FOS", "G0S2", "HBB", "HIST1H4C", "IER2", "IGKC", "IGLC2", "JUN", "JUNB", "KLRB1", "MT2A", "RPS26", "RPS4Y1", "TRBC1", "TUBA1B", "TUBB"))
  # 
  # excluded.genes <- c(excluded.genes,c("XIST", "MALAT1", "NEAT1"))
  # 
  # seurat.mc$innerNormVar <- mc_INV(pbmc[rowSums(pbmc[['RNA']]@counts)>40& ! rownames(pbmc[["RNA"]]) %in% excluded.genes,], membership.df)
  
  
  seurat.mc <- metacellCompactnessSeparation(sc.seurat = seurat.sc,SC_seurat = seurat.mc,dims = opt$RNAcomp)
  seurat.mc <- metacellCompactnessSeparation(sc.seurat = seurat.sc,SC_seurat = seurat.mc,preprocessing_method = "lsi", opt$ATACcomp)
  # Inner normalized variance
  seurat.mc$innerNormVar <- computeInnerNormVar(seurat = seurat.sc, memberships = seurat.mc@misc$membership)
  # Silhouette
  seurat.mc <- metacellSilhouetteSeurat(seurat.mc = seurat.mc,seurat = seurat.sc,reduction.name.sc = "pca",single.cells.res = T,dims =  opt$RNAcomp)
  seurat.mc <- metacellSilhouetteSeurat(seurat.mc = seurat.mc,seurat = seurat.sc,reduction.name.sc = "lsi",single.cells.res = T, dims = opt$ATACcomp) 
  
  
}



## Cell cycle analysis
DefaultAssay(seurat.mc) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat.mc <- CellCycleScoring(seurat.mc, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

if (!opt$useSize) {
  seurat.mc$size <- 1
}
pdf(paste0(opt$outdir,"/cellCyclePLot.pdf"))
plot(ggplot2::ggplot(seurat.mc@meta.data,aes(x=size,y=G2M.Score,color = Phase)) + geom_point())
dev.off()


if (!is.null(opt$pwmChromVAR)) {
  ## TF expression / TF motif accessibility correlation analyzis
  DefaultAssay(seurat.mc) <- peakAssay
  detectedMotifs <- rownames(seurat.mc@assays$chromvar) #motif
  relatedTFs <- ConvertMotifID(seurat.mc, id = detectedMotifs) #tf

  feature_x <- names(rowSums(seurat.mc@assays$RNA[rownames(seurat.mc@assays$RNA) %in% relatedTFs,]) > 0)[rowSums(seurat.mc@assays$RNA[rownames(seurat.mc@assays$RNA) %in% relatedTFs,]) > 0] # expressed tfs

  feature_y <- ConvertMotifID(seurat.mc, name = feature_x) # related motifs

  DefaultAssay(seurat.mc) <- "RNA"


  # DefaultAssay(seurat.mc.multi) <- "ACTIVITY"
  # seurat.mc.multi <- NormalizeData(seurat.mc.multi)

  chromVarCorrTable <- supercell_FeatureFeaturePlot_Seurat(seurat.mc = seurat.mc,
                                                           is.normalized = T,
                                                           cluster = "orig.ident",
                                                           assays = c("RNA","chromvar"),
                                                           method = "spearman",
                                                           feature_x = feature_x,
                                                           feature_y = feature_y,
                                                           use.size = opt$useSize,
                                                           plot = F)

  chromVarCorrTablePearson <- supercell_FeatureFeaturePlot_Seurat(seurat.mc = seurat.mc,
                                                                  is.normalized = T,
                                                                  cluster = "orig.ident",
                                                                  assays = c("RNA","chromvar"),
                                                                  method = "pearson",
                                                                  feature_x = feature_x,
                                                                  feature_y = feature_y,
                                                                  use.size = opt$useSize,
                                                                  plot = F)

  write.csv(chromVarCorrTable,paste0(opt$outdir,"/chromVarCorrTable.csv"))
  write.csv(chromVarCorrTablePearson,paste0(opt$outdir,"/chromVarCorrTablePearson.csv"))
}




saveRDS(seurat.mc,paste0(opt$outdir,"/seurat.multiome.mc.activities.rds"))

write.csv(corrTablePearson,paste0(opt$outdir,"/corrTablePearson.csv"))
write.csv(corrTable,paste0(opt$outdir,"/corrTable.csv"))

