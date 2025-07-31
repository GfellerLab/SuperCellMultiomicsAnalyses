# Load Packages
library(Pando)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)
library(purrr)
library("stringr")
library(parallel)
library(getopt)
library(dplyr)
library(SuperCell)
library(Signac)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'mcSeurat',  'm', 1, "character", "path to metacell data",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'ncores', "w", 1, "numeric", "number of core to use",
  "singleCellSeurat", "s", 1, "character", "Path to singlecell data",
  'mc.method',  'x', 1, "character", "metacell pando method",
  'membership',  'b', 1, "character", "metacell membership"
), byrow=TRUE, ncol=5)

opt = getopt(spec)

registerDoParallel(opt$ncores)

source("R/functions/pando_functions.r")


# Load sc or mc data ------------------------------------------------------


if(opt$mcSeurat != opt$singleCellSeurat){

  # Load single-cell object to get the same variable genes as the one used for single-cell data
  sc.seurat <- readRDS(opt$singleCellSeurat)
  sc.seurat <- DietSeurat(object = sc.seurat, assays = c("RNA", "ATAC"))
  gc()

  # Define most variable genes (used as targets later)
  DefaultAssay(sc.seurat) <- "RNA"
  sc.seurat <- NormalizeData(sc.seurat, assay='RNA')
  sc.seurat <- FindVariableFeatures(sc.seurat, assay='RNA', nfeatures = 4000)
  hvgs <- Seurat::VariableFeatures(sc.seurat)

  mc.use = T
  if(opt$mc.method == "aggregate"){

    # load membership
    membership_df <- read.csv(opt$membership, header = T, row.names = 1)
    membership <- membership_df$seurat.mc.multi.membership
    names(membership) <- rownames(membership_df)

    # add membership to the metadata
    seurat_object <- sc.seurat
    seurat_object$mc_membership <- membership[colnames(seurat_object)]
    rm(sc.seurat)
    gc()


  }else if(opt$mc.method %in% c("LSI1_norm", "LSI1_norm_svy")){
    rm(sc.seurat)
    gc()

    seurat_object <- readRDS(opt$mcSeurat)
    DefaultAssay(seurat_object) <- "RNA"
    seurat_object <- NormalizeData(seurat_object, assay='RNA')

    DefaultAssay(seurat_object) <- "ATAC"
    seurat_object <- Signac::RunTFIDF(seurat_object, assay='ATAC', method = 1)

  }else if(opt$mc.method %in% c("LSI2_norm", "LSI2_norm_svy")){
    rm(sc.seurat)
    gc()

    seurat_object <- readRDS(opt$mcSeurat)
    DefaultAssay(seurat_object) <- "RNA"
    seurat_object <- NormalizeData(seurat_object, assay='RNA')

    DefaultAssay(seurat_object) <- "ATAC"
    seurat_object <- Signac::RunTFIDF(seurat_object, assay='ATAC', method = 2)

  }else if(opt$mc.method  %in% c("meanSC_norm", "meanSC_norm_svy")){


    membership_df <- read.csv(opt$membership, header = T, row.names = 1)
    membership <- membership_df$seurat.mc.multi.membership
    names(membership) <- rownames(membership_df)
    seurat_object <- SCimplify_for_Seurat_v5(seurat = sc.seurat,
                                             membership = membership[colnames(sc.seurat)] ,
                                             assay = c("RNA", "ATAC"),
                                             avg.in.data = T,
                                             return.seurat = T)
    rm(sc.seurat)
    gc()

  }

}else{
  mc.use = F
  seurat_object <- readRDS(opt$singleCellSeurat)
  seurat_object <- DietSeurat(object = seurat_object, assays = c("RNA", "ATAC"))
  gc()

  # Define most variable genes (used as targets later)
  DefaultAssay(seurat_object) <- "RNA"
  seurat_object <- NormalizeData(seurat_object, assay='RNA')
  seurat_object <- FindVariableFeatures(seurat_object, assay='RNA', nfeatures = 4000)
  hvgs <- VariableFeatures(seurat_object)

}


# Get motif data ----------------------------------------------------------
data('phastConsElements20Mammals.UCSC.hg38')

# Initiate GRN object and select candidate regions ------------------------
grn_object <- initiate_grn(seurat_object,
                           peak_assay = "ATAC",
                           rna_assay = "RNA",
                           regions = phastConsElements20Mammals.UCSC.hg38)

data('motifs')
data('motif2tf')

motif2tf_use <- motif2tf %>%
  filter(tf %in%  hvgs)
motifs_use <- motifs[unique(motif2tf_use$motif)]
motif2tf_use

grn_object <- find_motifs(
  grn_object,
  pfm = motifs_use,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  motif_tfs = motif2tf_use
)


# Infer gene regulatory network (regions near genes are selected) ---------

methods = c("sc", "aggregate", "LSI1_norm", "LSI1_norm_svy", "LSI2_norm", "LSI2_norm_svy", "meanSC_norm", "meanSC_norm_svy")
model = c("glm","glm", "glm", "svyglm", "glm", "svyglm", "glm", "svyglm")
names(model) <- methods


if(mc.use){
  if(opt$mc.method == "aggregate"){
    grn_object <- infer_grn.GRNData(
      grn_object,
      aggregate_peaks_col = "mc_membership",
      aggregate_rna_col = "mc_membership",
      weights = NULL,
      peak_to_gene_method = 'GREAT',
      genes = hvgs,
      parallel = T
    )

  }else{
    grn_object <- infer_grn.GRNData(
      object = grn_object,
      method = model[opt$mc.method],
      weights = grn_object@data$size,
      peak_to_gene_method = 'GREAT',
      genes = hvgs,
      parallel = T
    )
  }

}else{
  grn_object <- infer_grn(grn_object,
                          genes = hvgs,
                          peak_to_gene_method = 'GREAT',
                          method = model[opt$mc.method], parallel = T)
}

grn_object@data <- DietSeurat(grn_object@data, assays = c("RNA", "ATAC"))
saveRDS(grn_object, file = paste0(opt$outdir, "grn_object_",opt$mc.method,".rds"))
