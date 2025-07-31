library(reticulate)
use_python("/opt/conda/envs/MetacellAnalysisToolkit/bin/python", required = TRUE)
library(Seurat)
library(dplyr)
library(getopt)
library(dplyr)
library(SuperCell)
library(ggplot2)


spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds) or dataset name (from SeuratData package)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "RNAcomp", "p", 1, "character", "range of RNA components to consider for metacell identification  (eg 1:50 for RNA pca)",
  "ADTcomp", "q", 1, "character", "range of ADT components to consider for metacell identification (eg 2:50 for ADT pca)",
  "qcRNAcomp", "v", 1, "character", "range of RNA components to consider for metacell checking  (eg 1:50 for RNA pca)",
  "qcADTcomp", "w", 1, "character", "range of ADT components to consider for metacell checking (eg 2:50 for ADT pca)",
  "RNAassay", "r", 1, "character", "RNA assay name with computed pca (default RNA)",
  "ADTassay", "a", 1, "character", "ADT assay name with computed apca (default ATAC)" ,
  "k.wnn", "k", 1, "numeric", "k for the knn used in the wnn used for metacell identification",
  "gamma", "g", 1, "numeric", "gamma used for metacell identificaiton",
  "kernel", "e", 1, "logical", "whether to use a kernel",
  "randomMetacells", "d", 1, "logical","construct random metacell at the specified gamma",
  "memberships", "c", 1 , "character", "memberships already computed (e.g. with SEACells) csv table: 1st column cell name, 2nd column metacell name",
  "prefixMembership", "x", 1, "character", "metacell name prefix memebership to disccard (e.g. SEACell-)",
  "inputSeuratMetacell", "m", 1, "character", "seurat metacell object if available to rescale directly",
  "returnRes", "s", 1, "character", "how to save the SuperCell results (seurat or SuperCellHierarchy, no save by default)",
  "doAnalysis", "y", 1, "character", "Set a method name (SuperCell_Multi) will make inter modality correlation and compactness/separation analysis)",
  "geneProtein", "z", 1, "character", "gene protein correspondance files for inter modalities correlation",
  "pythonSeacellEnv", "t", 1, "character", "python path for seacell env to compute compactness/separation"

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
# opt$RNAcomp <- "1:50"
# opt$ADTcomp <- "A:50"
# opt$minCutOff <- "q0"
# opt$RNAnormalization <- "SCTransform"

print(opt)


if (is.null(opt$RNAassay)) {
  opt$RNAassay <- "RNA"
}

if (is.null(opt$ADTassay)) {
  opt$ADTassay <- "ADT"
}

if (is.null(opt$randomMetacells)) {
  opt$randomMetacells <- F
}



if (!is.null(opt$RNAcomp)) {
  ci <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$RNAcomp,split = ":")[[1]][2])
  opt$RNAcomp <- c(ci:cf)
}



if (!is.null(opt$ADTcomp)) {
  ci <- as.numeric(strsplit(opt$ADTcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$ADTcomp,split = ":")[[1]][2])
  opt$ADTcomp <- c(ci:cf)
}


if (!is.null(opt$qcRNAcomp)) {
  ci <- as.numeric(strsplit(opt$qcRNAcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$qcRNAcomp,split = ":")[[1]][2])
  opt$qcRNAcomp <- c(ci:cf)
}

if (!is.null(opt$qcADTcomp)) {
  ci <- as.numeric(strsplit(opt$qcADTcomp,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$qcADTcomp,split = ":")[[1]][2])
  opt$qcADTcomp <- c(ci:cf)
}

if (is.null(opt$outdir)) {
  opt$outdir <- "./"
}

print(opt)


dir.create(opt$outdir,recursive = T,showWarnings = F)


if (endsWith(opt$inputSeurat,'h5ad')) {
  return.seurat <- F
  file.name <- strsplit(opt$inputSeurat,split = ".h5ad")[[1]][1]
  adata <- anndata::read_h5ad(opt$inputSeurat)
  if(!is.null(adata$raw)) {
    counts <- Matrix::t(adata$raw$X)
    rownames(counts) <- rownames(adata$raw$var)

  } else{
    counts <- Matrix::t(adata$X)
    rownames(counts) <- rownames(adata$var)
  }
  colnames(counts) <- adata$obs_names

  if (!grepl(x = file.name,pattern = "ADT")) {
    embeddings <- adata$obsm$X_pca
    rownames(embeddings) <- adata$obs_names
    seurat <- CreateSeuratObject(counts = counts,meta.data = adata$obs,assay=opt$RNAassay)
    seurat[["pca"]] <- CreateDimReducObject(embeddings = embeddings,key = "PCA_",assay = opt$RNAassay)
  } else {
    seurat <- CreateSeuratObject(counts = counts,meta.data = adata$obs,assay =opt$ADTassay)
    embeddings <- adata$obsm$X_apca
    rownames(embeddings) <- adata$obs_names
    seurat[["apca"]] <- CreateDimReducObject(embeddings = embeddings,key = "PC_",assay = opt$ADTassay)
  }
} else {
  seurat <- readRDS(opt$inputSeurat)
}



if (is.null(opt$inputSeuratMetacell)) {

  if (opt$randomMetacells) {
    print("constructing random metacell")
    randomMC <- 1:floor(ncol(seurat)/opt$gamma)
    randomMemberships <-sample(c(randomMC, sample(randomMC, ncol(seurat)-length(randomMC), replace=TRUE)))
    names(randomMemberships) <- colnames(seurat)
    seurat.mc.multi <- SCimplify_for_Seurat_v5(seurat = seurat,
                                               membership=randomMemberships)
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
    print(head(memberships))
    seurat.mc.multi <- SCimplify_for_Seurat_v5(seurat = seurat,
                                               membership=memberships,
                                               fragmentFiles = fragmentFiles)
  }

  if (!is.null(opt$RNAcomp)|!is.null(opt$ADTcomp)) {

    if (!is.null(opt$RNAcomp)&!is.null(opt$ADTcomp)) {

      seurat.mc.multi <- SCimplify_for_Seurat_v5(
        seurat,
        assay = c(opt$RNAassay,opt$ADTassay),
        k.knn = opt$k.wnn,
        reduction = list("pca", "apca"),
        dims = list(opt$RNAcomp, opt$ADTcomp),
        graph.name = "knn",
        kernel = opt$kernel,
        gamma = opt$gamma
      )} else {
        if (is.null(opt$RNAcomp)&!is.null(opt$ADTcomp)) {
          print("identifying metacells on ADT modality")
          seurat.mc.multi <- SCimplify_for_Seurat_v5(
            seurat,
            assay = c(opt$ADTassay),
            k.knn = opt$k.wnn,
            reduction = list("apca"),
            dims = list(opt$ADTcomp),
            fragmentFiles = fragmentFiles,
            tmpPath =paste0(opt$outdir,"tmp/"),
            outputDirMcFragment = paste0(getwd(),"/",outputDirMcFragment),
            kernel = opt$kernel,
            gamma = opt$gamma
          )
        }

        if (!is.null(opt$RNAcomp)&is.null(opt$ADTcomp)) {
          print("identifying metacells on RNA modality")
          seurat.mc.multi <- SCimplify_for_Seurat_v5(
            seurat,
            assay = c(opt$RNAassay),
            k.knn = opt$k.wnn,
            reduction = list("pca"),
            dims = list(opt$RNAcomp),
            kernel = opt$kernel,
            gamma = opt$gamma)
        }
      }
  }
} else {

  if (endsWith(opt$inputSeuratMetacell, "SuperCellHierarchy.rds")) { # SC hierarchy object as input
    seurat.mc <- seurat
    seurat.mc@misc$metacells_hierarchy<- readRDS(opt$inputSeuratMetacell)$metacells_hierarchy
  } else {
    seurat.mc <- readRDS(opt$inputSeuratMetacell)
  }
  seurat.mc.multi <- SCimplify_for_Seurat_v5(seurat = seurat, seurat.mc = seurat.mc, gamma = opt$gamma)
}

if (!is.null(opt$doAnalysis)) {

  #source("R/functions/innerNormalizedVariance.R")
  source("R/functions/metacellCompactnessSeparation.R")

  input <- opt$doAnalysis

  DefaultAssay(seurat.mc.multi) <- "RNA"
  seurat.mc.multi <- NormalizeData(seurat.mc.multi)

  DefaultAssay(seurat.mc.multi) <- "ADT"
  seurat.mc.multi <- NormalizeData(seurat.mc.multi, normalization.method = 'CLR', margin = 2)

  gene_protein_final <- read.csv(opt$geneProtein)

  w.cor <- FeatureFeaturePlot.SuperCell(seurat.obj = seurat.mc.multi,
                                               feature.x = gene_protein_final$gene.name,
                                               feature.y = gene_protein_final$X.protein,
                                               assays = c("RNA","ADT"),
                                               is.normalized = T,
                                               plot = F)

  w.cor.pearson <- FeatureFeaturePlot.SuperCell(seurat.obj = seurat.mc.multi,
                                                       feature.x = gene_protein_final$gene.name,
                                                       feature.y = gene_protein_final$X.protein,
                                                       assays = c("RNA","ADT"),
                                                       is.normalized = T,
                                                       plot = F)

  w.cor$input <- input
  w.cor$gamma <- opt$gamma
  w.cor$origIdent <- seurat.mc.multi$orig.ident[1]

  w.cor.pearson$input <- input
  w.cor.pearson$gamma <- opt$gamma
  w.cor.pearson$origIdent <- seurat.mc.multi$orig.ident[1]

  library(reticulate)
  use_python(opt$pythonSeacellEnv)
  library(MetacellAnalysisToolkit)
  # we discard metacell outliers with na.exclude
  membership_df <- data.frame("membership" = paste0("Metacell_",na.exclude(seurat.mc.multi@misc$membership)),
                              row.names = colnames(seurat)[!is.na(seurat.mc.multi@misc$membership)])

  seurat.mc.multi$compactness_pca <- mc_compactness(cell.membership = membership_df,
                                                    sc.obj = seurat,
                                                    sc.reduction = "pca_diffusion")
  seurat.mc.multi$compactness_apca <- mc_compactness(cell.membership = membership_df,
                                                     sc.obj = seurat,
                                                     sc.reduction = "apca_diffusion")
  seurat.mc.multi$separation_pca <- mc_separation(cell.membership = membership_df,
                                                   sc.obj = seurat,
                                                   sc.reduction = "pca_diffusion")
  seurat.mc.multi$separation_apca <- mc_separation(cell.membership = membership_df,
                                    sc.obj = seurat,
                                    sc.reduction = "apca_diffusion")

  # print(is(seurat.mc.multi))
  message("computing ASW")

  seurat.mc.multi$asw_apca_diffusion <- mc_ASW(seurat.mc = seurat.mc.multi,
                                              seurat.sc = seurat,
                                              cell.membership = membership_df,
                                              reduction.name.sc = "apca_diffusion")

  # print(is(seurat.mc.multi))
  # if(opt$gamma > 49) {
  #   seurat.mc.multi$innerNormVar <- mc_INV(cell.membership = membership_df,
  #                                          sc.obj = seurat)
  # }

  seurat.mc.multi$asw_pca_diffusion <- mc_ASW(seurat.mc  = seurat.mc.multi,
                                              seurat.sc = seurat,
                                              cell.membership = membership_df,
                                              reduction.name.sc = "pca_diffusion")


  message("all metrics computed")

  ## Inner normalized variance

  #seurat.mc.multi$innerNormVar <- computeInnerNormVar(seurat = seurat, memberships = seurat.mc.multi@misc$membership)


  # ## Cell cycle analysis To remove cell cycle phase is not available for all datasets
  #
  # DefaultAssay(seurat.mc.multi) <- "RNA"
  #
  # s.genes <- cc.genes$s.genes
  # g2m.genes <- cc.genes$g2m.genes
  # seurat.mc.multi <- CellCycleScoring(seurat.mc.multi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  #
  # pdf(paste0(opt$outdir,"/cellCyclePLot.pdf"))
  # plot(ggplot2::ggplot(seurat.mc.multi@meta.data,aes(x=size,y=G2M.Score,color = Phase)) + geom_point())
  # plot(ggplot2::ggplot(seurat.mc.multi@meta.data,aes(x=orig.ident,y=Phase_purity)) + geom_boxplot())
  # dev.off()


  write.csv(seurat.mc.multi@meta.data,paste0(opt$outdir,"/metaData.csv"))
  #write.csv(w.cor,paste0(opt$outdir,"/corrTableSpearman.csv"))
  write.csv(w.cor.pearson,paste0(opt$outdir,"/corrTablePearson.csv"))


}

if (! is.null(opt$returnRes)) {
  if (opt$returnRes == "seurat") {
    saveRDS(seurat.mc.multi, paste0(opt$outdir,"/seurat.cite.mc.rds"))
  }
  if (opt$returnRes == "SuperCellHierarchy") {
    saveRDS(seurat.mc.multi@misc, paste0(opt$outdir,"/SuperCellHierarchy.rds"))
  }
}
