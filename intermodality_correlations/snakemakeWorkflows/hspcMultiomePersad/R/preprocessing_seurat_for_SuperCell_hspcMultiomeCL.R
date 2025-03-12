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
  'fragmentFile', "f",2,  "character", "Fragment file path (rep1+rep2)",
  "k.wnn", "k", 1, "numeric", "k for the knn used in the wnn analysis",
  "RNAnormalization", "a", 1, "character", "normalisation method for RNA (logNormalize or SCTransform)", 
  "nVarGenes", "v", 1, "numeric", "number of variable genes",
  'minCutOff', "c", 1, "character", "ATAC features selection cut off (default q0)"
), byrow=TRUE, ncol=5)

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code
# test
# opt <- list()
# opt$fragmentFile <- "~/Documents/multiomicsMetacells/multiome_PBMC_data/fragments_files/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
# opt$inputSeurat <- "pbmcMultiome"
# opt$outdir <- "output/correlationAnalyzis/pbmcMultiome/singlecell_analysis"
# frag.file <- opt$fragmentFile 

# opt$minCutOff <- "q0"
# opt$RNAnormalization <- "SCTransform"

if(is.null(opt$RNAnormalization)) {
  opt$RNAnormalization <- "logNormalize"
}

if(!is.null(opt$fragmentFile)) {
  opt$fragmentFile <- strsplit(opt$fragmentFile,"\\+")[[1]]
}

if (is.null(opt$minCutOff)) {
  opt$minCutOff <- "q0"
}


if (is.null(opt$outdir)) {
  opt$outdir <- "./"
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

# frag_rep1 <- CreateFragmentObject(path = opt$fragmentFile[1],cells = rep1)
# frag_rep2 <- CreateFragmentObject(path = opt$fragmentFile[2],cells = rep2)


# Cells(frag_rep1) <- rep1.ori
# Cells(frag_rep2) <- rep2.ori
# frag.file <- list(frag_rep1,frag_rep2)
chrom_assay <- CreateChromatinAssay(
  counts = atac.count,
  sep = c(":", "-"),
  genome = 'hg38',
  # fragments = frag.file,
  #min.cells = 10,
  annotation = annotations
)
hspc[["ATAC"]] <- chrom_assay







if (opt$RNAnormalization == "SCTransform") {
  rnaAssay = "SCT"
  DefaultAssay(hspc) <- "RNA"
  hspc <- SCTransform(hspc, verbose = FALSE,conserve.memory = T) %>% RunPCA() 
} else {
  rnaAssay = "RNA"
  hspc <- NormalizeData(hspc, verbose = FALSE) %>% FindVariableFeatures(nFeature = opt$nVarGenes) %>% ScaleData() %>% RunPCA()
}


DefaultAssay(hspc) <- "ATAC"
hspc <- RunTFIDF(hspc)
hspc <- FindTopFeatures(hspc, min.cutoff = opt$minCutOff)
hspc <- RunSVD(hspc)




saveRDS(hspc, paste0(opt$outdir,"/seurat_multimodal.rds"))
