setwd("input/prostateCancer10xMultiome/")
BiocManager::install("DiffBind")
getwd()
{library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(EnsDb.Mmusculus.v79)
  library(Signac)
  library(S4Vectors)
  library(patchwork)
  set.seed(1234)}

#WT 
{inputdata.10x.WT.1 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/WT_1/filtered_feature_bc_matrix.h5")
  inputdata.10x.WT.2 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/WT_2/filtered_feature_bc_matrix.h5")
  atac_counts.WT.1 <- inputdata.10x.WT.1$Peaks
  atac_counts.WT.2 <- inputdata.10x.WT.2$Peaks

  #2W 
  inputdata.10x.2W.1 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2W_1//filtered_feature_bc_matrix.h5")
  inputdata.10x.2W.2 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2W_2//filtered_feature_bc_matrix.h5")
  atac_counts.2W.1 <- inputdata.10x.2W.1$Peaks
  atac_counts.2W.2 <- inputdata.10x.2W.2$Peaks

  #1M
  inputdata.10x.1M.1 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/1M_1/filtered_feature_bc_matrix.h5")
  inputdata.10x.1M.2 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/1M_2/filtered_feature_bc_matrix.h5")
  atac_counts.1M.1 <- inputdata.10x.1M.1$Peaks
  atac_counts.1M.2 <- inputdata.10x.1M.2$Peaks

  #2.5M 
  inputdata.10x.2.5M.1 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_1/filtered_feature_bc_matrix.h5")
  inputdata.10x.2.5M.2 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_2/filtered_feature_bc_matrix.h5")
  inputdata.10x.2.5M.3 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_3/filtered_feature_bc_matrix.h5")
  inputdata.10x.2.5M.4 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_4/filtered_feature_bc_matrix.h5")
  atac_counts.2.5M.1 <- inputdata.10x.2.5M.1$Peaks
  atac_counts.2.5M.2 <- inputdata.10x.2.5M.2$Peaks
  atac_counts.2.5M.3 <- inputdata.10x.2.5M.3$Peaks
  atac_counts.2.5M.4 <- inputdata.10x.2.5M.4$Peaks

  #3.5M 
  inputdata.10x.3.5M.1 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/3_5M_1/filtered_feature_bc_matrix.h5")
  inputdata.10x.3.5M.2 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/3_5M_2/filtered_feature_bc_matrix.h5")
  atac_counts.3.5M.1 <- inputdata.10x.3.5M.1$Peaks
  atac_counts.3.5M.2 <- inputdata.10x.3.5M.2$Peaks

  #4.5M 
  inputdata.10x.4.5M.1 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/4_5M_1/filtered_feature_bc_matrix.h5")
  inputdata.10x.4.5M.2 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/4_5M_2/filtered_feature_bc_matrix.h5")
  atac_counts.4.5M.1 <- inputdata.10x.4.5M.1$Peaks
  atac_counts.4.5M.2 <- inputdata.10x.4.5M.2$Peaks

  #RC
  inputdata.10x.RC.1 <- Read10X_h5("sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/6M/filtered_feature_bc_matrix.h5")
  atac_counts.RC.1 <- inputdata.10x.RC.1$Peaks}

#######Define frag.files
{frag.file.WT.1 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/WT_1/atac_fragments.tsv.gz"
  frag.file.WT.2 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/WT_2/atac_fragments.tsv.gz"
  frag.file.2W.1 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2W_1/atac_fragments.tsv.gz"
  frag.file.2W.2 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2W_2/atac_fragments.tsv.gz"
  frag.file.1M.1 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/1M_1/atac_fragments.tsv.gz"
  frag.file.1M.2 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/1M_2/atac_fragments.tsv.gz"
  frag.file.2.5M.1 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_1/atac_fragments.tsv.gz"
  frag.file.2.5M.2 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_2/atac_fragments.tsv.gz"
  frag.file.2.5M.3 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_3/atac_fragments.tsv.gz"
  frag.file.2.5M.4 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/2_5M_4/atac_fragments.tsv.gz"
  frag.file.3.5M.1 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/3_5M_1/atac_fragments.tsv.gz"
  frag.file.3.5M.2 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/3_5M_2/atac_fragments.tsv.gz"
  frag.file.4.5M.1 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/4_5M_1/atac_fragments.tsv.gz"
  frag.file.4.5M.2 <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/4_5M_2/atac_fragments.tsv.gz"
  frag.file.RC <- "sibcb2/gaodonglab2/lifei/NEPC/sc_REG/sc-RNA-ATAC/cell_ranger_outs/6M/atac_fragments.tsv.gz"}
#Create Seurat object
#construct atac seurat object
#atac.WT.1
{{grange.counts.WT.1 <- StringToGRanges(rownames(atac_counts.WT.1), sep = c(":", "-"))
grange.use.WT.1 <- seqnames(grange.counts.WT.1) %in% standardChromosomes(grange.counts.WT.1)
atac_counts.WT.1 <- atac_counts.WT.1[as.vector(grange.use.WT.1), ]
annotations.WT.1 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations.WT.1) <- 'UCSC'
genome(annotations.WT.1) <- "mm10"

chrom_assay.WT.1 <- CreateChromatinAssay(
  counts = atac_counts.WT.1,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file.WT.1,
  min.cells = 10,
  annotation = annotations.WT.1
)}
  #atac.WT.2
  {grange.counts.WT.2 <- StringToGRanges(rownames(atac_counts.WT.2), sep = c(":", "-"))
    grange.use.WT.2 <- seqnames(grange.counts.WT.2) %in% standardChromosomes(grange.counts.WT.2)
    atac_counts.WT.2 <- atac_counts.WT.2[as.vector(grange.use.WT.2), ]
    annotations.WT.2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.WT.2) <- 'UCSC'
    genome(annotations.WT.2) <- "mm10"
    
    chrom_assay.WT.2 <- CreateChromatinAssay(
      counts = atac_counts.WT.2,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.WT.2,
      min.cells = 10,
      annotation = annotations.WT.2
    )}
  
  #atac.2W.1
  {grange.counts.2W.1 <- StringToGRanges(rownames(atac_counts.2W.1), sep = c(":", "-"))
    grange.use.2W.1 <- seqnames(grange.counts.2W.1) %in% standardChromosomes(grange.counts.2W.1)
    atac_counts.2W.1 <- atac_counts.2W.1[as.vector(grange.use.2W.1), ]
    annotations.2W.1 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.2W.1) <- 'UCSC'
    genome(annotations.2W.1) <- "mm10"
    
    chrom_assay.2W.1 <- CreateChromatinAssay(
      counts = atac_counts.2W.1,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.2W.1,
      min.cells = 10,
      annotation = annotations.2W.1
    )}
  
  
  #atac.2W.2
  {grange.counts.2W.2 <- StringToGRanges(rownames(atac_counts.2W.2), sep = c(":", "-"))
    grange.use.2W.2 <- seqnames(grange.counts.2W.2) %in% standardChromosomes(grange.counts.2W.2)
    atac_counts.2W.2 <- atac_counts.2W.2[as.vector(grange.use.2W.2), ]
    annotations.2W.2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.2W.2) <- 'UCSC'
    genome(annotations.2W.2) <- "mm10"
    
    chrom_assay.2W.2 <- CreateChromatinAssay(
      counts = atac_counts.2W.2,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.2W.2,
      min.cells = 10,
      annotation = annotations.2W.2
    )}
  
  #atac.1M.1
  {grange.counts.1M.1 <- StringToGRanges(rownames(atac_counts.1M.1), sep = c(":", "-"))
    grange.use.1M.1 <- seqnames(grange.counts.1M.1) %in% standardChromosomes(grange.counts.1M.1)
    atac_counts.1M.1 <- atac_counts.1M.1[as.vector(grange.use.1M.1), ]
    annotations.1M.1 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.1M.1) <- 'UCSC'
    genome(annotations.1M.1) <- "mm10"
    
    chrom_assay.1M.1 <- CreateChromatinAssay(
      counts = atac_counts.1M.1,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.1M.1,
      min.cells = 10,
      annotation = annotations.1M.1
    )}
  
  #atac.1M.2
  {grange.counts.1M.2 <- StringToGRanges(rownames(atac_counts.1M.2), sep = c(":", "-"))
    grange.use.1M.2 <- seqnames(grange.counts.1M.2) %in% standardChromosomes(grange.counts.1M.2)
    atac_counts.1M.2 <- atac_counts.1M.2[as.vector(grange.use.1M.2), ]
    annotations.1M.2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.1M.2) <- 'UCSC'
    genome(annotations.1M.2) <- "mm10"
    
    chrom_assay.1M.2 <- CreateChromatinAssay(
      counts = atac_counts.1M.2,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.1M.2,
      min.cells = 10,
      annotation = annotations.1M.2
    )}
  #atac.2.5M.1
  {grange.counts.2.5M.1 <- StringToGRanges(rownames(atac_counts.2.5M.1), sep = c(":", "-"))
    grange.use.2.5M.1 <- seqnames(grange.counts.2.5M.1) %in% standardChromosomes(grange.counts.2.5M.1)
    atac_counts.2.5M.1 <- atac_counts.2.5M.1[as.vector(grange.use.2.5M.1), ]
    annotations.2.5M.1 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.2.5M.1) <- 'UCSC'
    genome(annotations.2.5M.1) <- "mm10"
    
    chrom_assay.2.5M.1 <- CreateChromatinAssay(
      counts = atac_counts.2.5M.1,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.2.5M.1,
      min.cells = 10,
      annotation = annotations.2.5M.1
    )}
  
  #atac.2.5M.2
  {grange.counts.2.5M.2 <- StringToGRanges(rownames(atac_counts.2.5M.2), sep = c(":", "-"))
    grange.use.2.5M.2 <- seqnames(grange.counts.2.5M.2) %in% standardChromosomes(grange.counts.2.5M.2)
    atac_counts.2.5M.2 <- atac_counts.2.5M.2[as.vector(grange.use.2.5M.2), ]
    annotations.2.5M.2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.2.5M.2) <- 'UCSC'
    genome(annotations.2.5M.2) <- "mm10"
    
    chrom_assay.2.5M.2 <- CreateChromatinAssay(
      counts = atac_counts.2.5M.2,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.2.5M.2,
      min.cells = 10,
      annotation = annotations.2.5M.2
    )}
  
  #atac.2.5M.3
  {grange.counts.2.5M.3 <- StringToGRanges(rownames(atac_counts.2.5M.3), sep = c(":", "-"))
    grange.use.2.5M.3 <- seqnames(grange.counts.2.5M.3) %in% standardChromosomes(grange.counts.2.5M.3)
    atac_counts.2.5M.3 <- atac_counts.2.5M.3[as.vector(grange.use.2.5M.3), ]
    annotations.2.5M.3 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.2.5M.3) <- 'UCSC'
    genome(annotations.2.5M.3) <- "mm10"
    
    chrom_assay.2.5M.3 <- CreateChromatinAssay(
      counts = atac_counts.2.5M.3,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.2.5M.3,
      min.cells = 10,
      annotation = annotations.2.5M.3
    )}
  #atac.2.5M.4
  {grange.counts.2.5M.4 <- StringToGRanges(rownames(atac_counts.2.5M.4), sep = c(":", "-"))
    grange.use.2.5M.4 <- seqnames(grange.counts.2.5M.4) %in% standardChromosomes(grange.counts.2.5M.4)
    atac_counts.2.5M.4 <- atac_counts.2.5M.4[as.vector(grange.use.2.5M.4), ]
    annotations.2.5M.4 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.2.5M.4) <- 'UCSC'
    genome(annotations.2.5M.4) <- "mm10"
    
    chrom_assay.2.5M.4 <- CreateChromatinAssay(
      counts = atac_counts.2.5M.4,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.2.5M.4,
      min.cells = 10,
      annotation = annotations.2.5M.4
    )}
  
  #atac.3.5M.1
  {grange.counts.3.5M.1 <- StringToGRanges(rownames(atac_counts.3.5M.1), sep = c(":", "-"))
    grange.use.3.5M.1 <- seqnames(grange.counts.3.5M.1) %in% standardChromosomes(grange.counts.3.5M.1)
    atac_counts.3.5M.1 <- atac_counts.3.5M.1[as.vector(grange.use.3.5M.1), ]
    annotations.3.5M.1 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.3.5M.1) <- 'UCSC'
    genome(annotations.3.5M.1) <- "mm10"
    
    chrom_assay.3.5M.1 <- CreateChromatinAssay(
      counts = atac_counts.3.5M.1,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.3.5M.1,
      min.cells = 10,
      annotation = annotations.3.5M.1
    )}
  
  #atac.3.5M.2
  {grange.counts.3.5M.2 <- StringToGRanges(rownames(atac_counts.3.5M.2), sep = c(":", "-"))
    grange.use.3.5M.2 <- seqnames(grange.counts.3.5M.2) %in% standardChromosomes(grange.counts.3.5M.2)
    atac_counts.3.5M.2 <- atac_counts.3.5M.2[as.vector(grange.use.3.5M.2), ]
    annotations.3.5M.2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.3.5M.2) <- 'UCSC'
    genome(annotations.3.5M.2) <- "mm10"
    
    chrom_assay.3.5M.2 <- CreateChromatinAssay(
      counts = atac_counts.3.5M.2,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.3.5M.2,
      min.cells = 10,
      annotation = annotations.3.5M.2
    )}
  
  #atac.4.5M.1
  {grange.counts.4.5M.1 <- StringToGRanges(rownames(atac_counts.4.5M.1), sep = c(":", "-"))
    grange.use.4.5M.1 <- seqnames(grange.counts.4.5M.1) %in% standardChromosomes(grange.counts.4.5M.1)
    atac_counts.4.5M.1 <- atac_counts.4.5M.1[as.vector(grange.use.4.5M.1), ]
    annotations.4.5M.1 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.4.5M.1) <- 'UCSC'
    genome(annotations.4.5M.1) <- "mm10"
    
    chrom_assay.4.5M.1 <- CreateChromatinAssay(
      counts = atac_counts.4.5M.1,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.4.5M.1,
      min.cells = 10,
      annotation = annotations.4.5M.1
    )}
  
  #atac.4.5M.2
  {grange.counts.4.5M.2 <- StringToGRanges(rownames(atac_counts.4.5M.2), sep = c(":", "-"))
    grange.use.4.5M.2 <- seqnames(grange.counts.4.5M.2) %in% standardChromosomes(grange.counts.4.5M.2)
    atac_counts.4.5M.2 <- atac_counts.4.5M.2[as.vector(grange.use.4.5M.2), ]
    annotations.4.5M.2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.4.5M.2) <- 'UCSC'
    genome(annotations.4.5M.2) <- "mm10"
    
    chrom_assay.4.5M.2 <- CreateChromatinAssay(
      counts = atac_counts.4.5M.2,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.4.5M.2,
      min.cells = 10,
      annotation = annotations.4.5M.2
    )}
  
  #atac.RC
  {grange.counts.RC <- StringToGRanges(rownames(atac_counts.RC.1), sep = c(":", "-"))
    grange.use.RC <- seqnames(grange.counts.RC) %in% standardChromosomes(grange.counts.RC)
    atac_counts.RC.1 <- atac_counts.RC.1[as.vector(grange.use.RC), ]
    annotations.RC <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.RC) <- 'UCSC'
    genome(annotations.RC) <- "mm10"
    
    chrom_assay.RC <- CreateChromatinAssay(
      counts = atac_counts.RC.1,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.RC,
      min.cells = 10,
      annotation = annotations.RC
    )}}
####release storage 
rm(atac_counts.1M.1,atac_counts.1M.2,atac_counts.2.5M.1,atac_counts.2.5M.2,atac_counts.2.5M.3,atac_counts.2.5M.4)
rm(atac_counts.2W.1,atac_counts.2W.2,atac_counts.3.5M.1,atac_counts.3.5M.2,atac_counts.4.5M.1,atac_counts.4.5M.2)
rm(atac_counts.WT.1,atac_counts.WT.2,atac_counts.RC.1)
gc()

##############construct consensus_peaksets
#BiocManager::install("DiffBind")
#library(amap)
library(DiffBind)
data<-"SampleSheet.csv"
dbObj <- dba(sampleSheet=data)
consensus_peakset <- dba.peakset(dbObj,dbObj$masks$Consensus,bRetrieve=TRUE)
dbObj_consensus1 <- dba.peakset(dbObj, consensus = DBA_CONDITION, minOverlap=0.5)#
dim(dbObj_consensus1)
consensus_peakset <- dba.peakset(dbObj_consensus1,dbObj_consensus1$masks$Consensus,bRetrieve=TRUE)
consensus_peakset<-consensus_peakset[1:105020,]
length(consensus_peakset)
class(consensus_peakset)
write.csv(consensus_peakset,"peakset_minOverlap_0.5_remove.csv")


















