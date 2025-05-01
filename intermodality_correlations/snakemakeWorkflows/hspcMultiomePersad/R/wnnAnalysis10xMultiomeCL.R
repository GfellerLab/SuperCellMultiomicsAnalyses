library(reticulate)
use_python("/opt/conda/envs/MetacellAnalysisToolkit/bin/python", required = TRUE)
library(SeuratData)
library(Seurat)
library(Signac)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(getopt)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputRNA',  'i', 1, "character", "REQUIRED : rna adata object (.h5ad)",
  'inputATAC',  'j', 1, "character", "REQUIRED : atac adata object (.h5ad)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'fragmentFile', "f",2,  "character", "Fragment file path",
  "RNAcomp", "p", 1, "character", "range of components to consider for wnn analysis (eg 1:50 for RNA pca)",
  "ATACcomp", "q", 1, "character", "range of components to consider for wnn analysis (eg 2:50 for ATAC pca)",
  "k.wnn", "k", 1, "numeric", "k for the knn used in the wnn analysis",
  "RNAnormalization", "a", 1, "character", "normalisation method for RNA (logNormalize or SCTransform)",
  "nVarGenes", "v", 1, "numeric", "number of variable genes",
  'minCutOff', "c", 1, "character", "ATAC features selection cut off (default q0)",
  "threads", "t", 1, "numeric", "number of threads to merge fragments files"
), byrow=TRUE, ncol=5)

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code
# test
# machine_path <- "/mnt/curnagl/"
# machine_path <- ""
# opt <- list()
# opt$fragmentFile <- paste0(machine_path, "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/SuperCellMultiomicsAnalyses/input/hspcMultiomePersad/BM_CD34_Rep1_atac_fragments.tsv.gz+/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/SuperCellMultiomicsAnalyses/input/hspcMultiomePersad/BM_CD34_Rep2_atac_fragments.tsv.gz")
# opt$inputRNA <- paste0(machine_path, "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/SuperCellMultiomicsAnalyses/input/hspcMultiomePersad/cd34_multiome_rna.h5ad")
# opt$inputATAC <- paste0(machine_path, "/work/FAC/FBM/LLB/dgfeller/scrnaseq/agabrie4/supercellV2/SuperCellMultiomicsAnalyses/input/hspcMultiomePersad/cd34_multiome_atac.h5ad")
# # frag.file <- opt$fragmentFile
# opt$RNAcomp <- "1:50"
# opt$ATACcomp <- "2:50"
# opt$minCutOff <- "q0"
# opt$RNAnormalization <- "SCTransform"
# opt$threads <- 8

if(is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "logNormalize"
}

if(!is.null(opt$fragmentFile)) {
  opt$fragmentFile <- strsplit(opt$fragmentFile,"\\+")[[1]]
}



if (is.null(opt$minCutOff)) {
  opt$minCutOff <- "q0"
}

if (is.null(opt$RNAcomp)) {
  opt$RNAcomp <- c(1:30)
} else {
  ci <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][2])
  opt$RNAcomp <- c(ci:cf)
}

if (is.null(opt$ATACcomp)) {
  opt$ATACcomp <- c(1:30)
} else {
  ci <- as.numeric(strsplit(opt$ATACcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$ATACcomp,split = ":")[[1]][2])
  opt$ATACcomp <- c(ci:cf)
}

if (is.null(opt$outdir)) {
  opt$outdir <- "./"
}

if(is.null(opt$threads)){
  opt$threads = parallel::detectCores()-4
}

dir.create(opt$outdir,recursive = T,showWarnings = F)




rna <- anndata::read_h5ad(opt$inputRNA)
atac <- anndata::read_h5ad(opt$inputATAC)

rna.count <- Matrix::t(rna$raw$X)
rownames(rna.count) <- rna$var_names
colnames(rna.count) <- rna$obs_names

hspc <- CreateSeuratObject(counts = rna.count,meta.data = rna$obs)

rep1.ori <- colnames(hspc)[grepl(x = colnames(hspc),pattern = "rep1")]
rep1 <- sub(x= rep1.ori,"cd34_multiome_rep1#",replacement = "")

rep2.ori <- colnames(hspc)[grepl(x = colnames(hspc),pattern = "rep2")]
rep2 <- sub(x= rep2.ori,"cd34_multiome_rep2#",replacement = "")
names(rep2) <- rep2.ori
names(rep1) <- rep1.ori

table(rep1 %in% rep2)

atac.count <- Matrix::t(atac$X)

grange.counts <- StringToGRanges(rownames(atac.count), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac.count <- atac.count[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

fragments.list <- list()
fragments.list[["rep1"]] <- CreateFragmentObject(path = opt$fragmentFile[1],cells = rep1)
fragments.list[["rep2"]] <- CreateFragmentObject(path = opt$fragmentFile[2],cells = rep2)

merge.fragmentsFiles <- function(fragments, output_path, output_prefix, threads){
  for(f in fragments){
    file = f@path
    # print(file)
    prefix <- unlist(strsplit(names(f@cells)[1],"#"))[1]
    system(paste0("gzip -dc ", file, "| awk 'BEGIN {FS=OFS=\"\\t\"} {print $1,$2,$3,", "\"", prefix,"_\"$4,$5}' - > ", output_path, "/", prefix, ".tsv"))
  }

  tsv_files <- list.files(path = output_path, pattern = ".tsv", full.names = T)
  tsv_files <- tsv_files[-which(grepl(".gz", tsv_files))]

  system(paste0("for file in ", paste(tsv_files, collapse = " "), "; do grep \"#\" \"$file\"; done > ", output_path,  "/header_lines.tsv"))
  system(paste0("for file in ", paste(tsv_files, collapse = " "), "; do grep -v \"#\" \"$file\"; done > ", output_path,  "/data_lines.tsv"))

  system(paste0("sort -k1,1V -k2,2n -k3,3n ", output_path,  "/data_lines.tsv", " > ", output_path,  "/sorted_data_lines.tsv"))
  system(paste0("cat ", output_path,  "/header_lines.tsv ", output_path,  "/sorted_data_lines.tsv > ", output_path,  "/", output_prefix, ".tsv"))

  system(paste0("bgzip -f -@ ", threads, " ", output_path,  "/", output_prefix, ".tsv"))

  system(paste0("tabix -p bed ", output_path,  "/", output_prefix, ".tsv.gz"))

  system(paste0("rm ", paste(c(tsv_files, paste0(output_path, "/", c("header_lines.tsv", "data_lines.tsv", "sorted_data_lines.tsv"))), collapse = " ") ))

  return( normalizePath(paste0(output_path,  "/", output_prefix, ".tsv.gz")) )
}

# Signac:::CallPeaks.Seurat()
frag.file <- merge.fragmentsFiles(fragments = fragments.list,
                                  output_path = dirname(opt$fragmentFile[1]),
                                  output_prefix = "combined_fragments",
                                  threads = opt$threads)

cell.names <- gsub("#", "_", colnames(hspc))
names(cell.names) <- colnames(hspc)
frag.obj <- CreateFragmentObject(path = frag.file, cells = cell.names)

# Cells(frag_rep1) <- rep1.ori
# Cells(frag_rep2) <- rep2.ori
# frag.file <- list(frag_rep1,frag_rep2)
chrom_assay <- CreateChromatinAssay(
  counts = atac.count,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.obj,
  #min.cells = 10,
  annotation = annotations
)
hspc[["ATAC"]] <- chrom_assay


# pbmc$coarse.annotation <- pbmc$seurat_annotations
#
# #pbmc$coarse.annotation[grepl(pattern = "B",x = pbmc$coarse.annotation)] <- "B"
#
# pbmc$coarse.annotation[grepl(pattern = "CD8 TEM",x = pbmc$coarse.annotation)] <- "CD8 Mem"
#
# pbmc$coarse.annotation[grepl(pattern = "CD4 TEM",x = pbmc$coarse.annotation)] <- "CD4 Mem"
# pbmc$coarse.annotation[grepl(pattern = "CD4 TCM",x = pbmc$coarse.annotation)] <- "CD4 Mem"




# #We define a color palette for this new annotations.
#
# color <- c("CD4 Naive"="#999999","NK"="#004949","CD8 Naive"="#009292","CD14 Mono"="#ff6db6",
#            "gdT"="#490092", "CD4 Mem"="#006ddb","cDC"="#b66dff","Treg"="#6db6ff",
#            "Intermediate B"="#b6dbff","Memory B"= "#8494FF","Naive B" = "#00A9FF",
#            "CD16 Mono"="#920000","HSPC"="#924900","CD8 Mem"="#db6d00","pDC"="#24ff24", "MAIT"="#ffff6d","Plasma"="#ffb6db")


## Analyzis of each modality separately

# Seurat scRNA-seq workflow

#As in Seurat tutorial for multimodal analyizis we use the SCTransform normalization for RNA data

if (opt$RNAnormalization == "SCTransform") {
  rnaAssay = "SCT"
  DefaultAssay(hspc) <- "RNA"
  hspc <- SCTransform(hspc, verbose = FALSE,conserve.memory = TRUE) %>% RunPCA()
} else {
  rnaAssay = "RNA"
  hspc <- NormalizeData(hspc, verbose = FALSE) %>% FindVariableFeatures(nFeature = opt$nVarGenes) %>% ScaleData() %>% RunPCA()
}


# Signac scATAC-seq workflow


# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
#grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
DefaultAssay(hspc) <- "ATAC"
hspc <- RunTFIDF(hspc)
hspc <- FindTopFeatures(hspc, min.cutoff = opt$minCutOff)
hspc <- RunSVD(hspc)




## Multimodal analyzis with Seurat

hspc <- RunUMAP(hspc,dims = opt$RNAcomp,
                reduction.name = 'umap.rna',
                reduction.key = 'rnaUMAP_')

hspc <- RunUMAP(hspc, reduction = 'lsi', dims = opt$ATACcomp, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


# UMAP results


p1 <- DimPlot(hspc, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("RNA")+ NoLegend()
p2 <- DimPlot(hspc, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("ATAC")+ NoLegend()
p3 <- DimPlot(hspc, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE, group.by = "sample") + ggplot2::ggtitle("RNA")+ NoLegend()
p4 <- DimPlot(hspc, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE, group.by = "sample") + ggplot2::ggtitle("ATAC")+ NoLegend()

pdf(file = paste0(opt$outdir,"/umap_single_modality.pdf"))
p1
p2
p3
p4
dev.off()

hspc <- FindMultiModalNeighbors(hspc, reduction.list = list("pca", "lsi"),k.nn = 30, dims.list = list(opt$RNAcomp, opt$ATACcomp))
hspc <- RunUMAP(hspc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",return.model = T)

#Comparison of UMAP results

p5 <- DimPlot(hspc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("WNN")
p6 <- DimPlot(hspc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE, group.by = "sample") + ggplot2::ggtitle("WNN")

pdf(file = paste0(opt$outdir,"/umap_wnn_analysis.pdf"))
p5
p6
dev.off()

hspc <- FindClusters(hspc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

pdf(file = paste0(opt$outdir,"/umap_wnn_analysis_clusters.pdf"))
DimPlot(hspc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggplot2::ggtitle("WNN")
dev.off()

# compute gene activity
#
DefaultAssay(hspc) <- "RNA"
hspc  <- NormalizeData(hspc)

genes <- rownames(hspc[["RNA"]]@layers$counts)[Matrix::rowSums(hspc[["RNA"]]@layers$counts)> ncol(hspc)*0.005]
write.table(genes,paste0(opt$outdir,"/selectedGenes.txt"))
#
DefaultAssay(hspc) <- "ATAC"
gene.activities <- GeneActivity(hspc,features =genes)


# add gene activities as a new assay
hspc[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

saveRDS(hspc, paste0(opt$outdir,"/seuratWNN.rds"))
