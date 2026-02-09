library(reticulate)
use_python("/opt/conda/envs/MetacellAnalysisToolkit/bin/python", required = TRUE)
library(SeuratData)
library(Seurat)
library(Signac)
library(dplyr)
# library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SuperCell)
library(getopt)
library(doParallel)
set.seed(2026)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED : seurat object (.rds) or dataset name (from SeuratData package)",
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
  "aggregateFragmentfile", "f", 0, "logical", 'whether to aggregate fragment file or not (default FALSE)',
  "returnMemberships", "b", 0, "logical", "wether to return only memberships (default seurat object with memberships inn misc slot)"

), byrow=TRUE, ncol=5)

opt = getopt(spec)


if(is.null(opt$returnMembership)) {
  return.seurat <- T
} else {
  return.seurat <- F
}

if (is.null(opt$RNAassay)) {
  opt$RNAassay <- "RNA"
}

if (is.null(opt$ATACassay)) {
  opt$ATACassay <- "ATAC"
}

if(is.null(opt$k.wnn)) {
  opt$k.wnn <- 30
}

if (is.null(opt$aggregateFragmentfile)) {
  opt$aggregateFragmentfile <- F
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
  if (!grepl(x = file.name,pattern = "ATAC")) {
    embeddings <- adata$obsm$X_pca
    rownames(embeddings) <- adata$obs_names
    seurat <- CreateSeuratObject(counts = counts,meta.data = adata$obs,assay=opt$RNAassay)
    seurat[["pca"]] <- CreateDimReducObject(embeddings = embeddings,key = "PCA_",assay = opt$RNAassay)
  } else {
    seurat <- CreateSeuratObject(counts = counts,meta.data = adata$obs,assay =opt$ATACassay)
    embeddings <- adata$obsm$X_lsi
    rownames(embeddings) <- adata$obs_names
    seurat[["lsi"]] <- CreateDimReducObject(embeddings = embeddings,key = "LSI_",assay = opt$ATACassay)
  }
} else {
  seurat <- readRDS(opt$inputSeurat)
}


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
    randomMemberships.v5 <- paste0("Metacell_",randomMemberships)

    if(opt$aggregateFragmentfile) {
      cell.names.fragments <- Fragments(seurat)[[1]]@cells
      names(randomMemberships.v5) <- cell.names.fragments[names(randomMemberships)]
    } else {
      names(randomMemberships.v5) <- names(randomMemberships)
    }

    seurat.mc.multi <- SCimplify_for_Seurat_v5(seurat = seurat,
                                               membership=randomMemberships)
    # fragmentFiles = fragmentFiles,
    # tmpPath =paste0(opt$outdir,"tmp/"),
    # outputDirMcFragment = paste0(getwd(),"/",outputDirMcFragment))
    mc.frag.path <- AggregateFragmentFile(input_file = fragmentFiles[["ATAC"]],
                                          tmp_path =paste0(opt$outdir,"/tmp/"),
                                          output_path = outputDirMcFragment,
                                          membership = randomMemberships.v5)

    frag.mc <- CreateFragmentObject(mc.frag.path,
                                    cells = colnames(seurat.mc.multi))

    Fragments(seurat.mc.multi[["ATAC"]]) <- frag.mc

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

    #AggregateExpression in SCimplify_for_Seurat_v5 remove cells with NA (MetaCell outliers) before aggregating expression
    seurat.mc.multi <- SCimplify_for_Seurat_v5(seurat = seurat,
                                               membership=memberships
    )


    memberships.v5 <- paste0("Metacell_",memberships)

    if(opt$aggregateFragmentfile) {
      cell.names.fragments <- Fragments(seurat)[[1]]@cells
      names(memberships.v5) <- cell.names.fragments[names(memberships)]
    } else {
      names(memberships.v5) <- names(memberships)
    }

    mc.frag.path <- AggregateFragmentFile(input_file = fragmentFiles[["ATAC"]],
                                          tmp_path =paste0(opt$outdir,"/tmp/"),
                                          output_path = outputDirMcFragment,
                                          membership = memberships.v5)

    frag.mc <- CreateFragmentObject(mc.frag.path,
                                    cells = colnames(seurat.mc.multi))

    Fragments(seurat.mc.multi[["ATAC"]]) <- frag.mc

  }

  if (!is.null(opt$RNAcomp)|!is.null(opt$ATACcomp)) {

    if (!is.null(opt$RNAcomp)&!is.null(opt$ATACcomp)) {

      seurat.mc.multi <- SCimplify_for_Seurat_v5(
        seurat,
        assay = c(opt$RNAassay,opt$ATACassay),
        k.knn = opt$k.wnn,
        reduction = list("pca", "lsi"),
        dims = list(opt$RNAcomp, opt$ATACcomp),
        graph.name = "knn",
        kernel = opt$kernel,
        gamma = opt$gamma,
        return.seurat = return.seurat
      )

      if (return.seurat) {
        memberships.v5 <- paste0("Metacell_",memberships)
        names(memberships.v5) <- names(memberships)

        mc.frag.path <- AggregateFragmentFile(input_file = fragmentFiles[["ATAC"]],
                                              tmp_path =paste0(opt$outdir,"/tmp/"),
                                              output_path = outputDirMcFragment,
                                              membership = memberships.v5)

        frag.mc <- CreateFragmentObject(mc.frag.path,
                                        cells = colnames(seurat.mc.multi))

        Fragments(seurat.mc.multi[["ATAC"]]) <- frag.mc

      }





    } else {
      if (is.null(opt$RNAcomp)&!is.null(opt$ATACcomp)) {
        print("identifying metacells on ATAC modality")
        seurat.mc.multi <- SCimplify_for_Seurat_v5(
          seurat,
          assay = c(opt$ATACassay),
          k.knn = opt$k.wnn,
          reduction = list("lsi"),
          dims = list(opt$ATACcomp),
          kernel = opt$kernel,
          gamma = opt$gamma,
          return.seurat = return.seurat
        )
        if (return.seurat) {
          memberships.v5 <- paste0("Metacell_",memberships)
          names(memberships.v5) <- names(memberships)

          mc.frag.path <- AggregateFragmentFile(input_file = fragmentFiles[["ATAC"]],
                                                tmp_path =paste0(opt$outdir,"/tmp/"),
                                                output_path = outputDirMcFragment,
                                                membership = memberships.v5)

          frag.mc <- CreateFragmentObject(mc.frag.path,
                                          cells = colnames(seurat.mc.multi))

          Fragments(seurat.mc.multi[["ATAC"]]) <- frag.mc

        }
      }

      if (!is.null(opt$RNAcomp)&is.null(opt$ATACcomp)) {
        print("identifying metacells on RNA modality")
        seurat.mc.multi <- SCimplify_for_Seurat_v5(
          seurat,
          assay = c(opt$RNAassay),
          k.knn = opt$k.wnn,
          reduction = list("pca"),
          dims = list(opt$RNAcomp),
          kernel = opt$kernel,
          gamma = opt$gamma,
          return.seurat = return.seurat)

        if (return.seurat) {
          memberships.v5 <- paste0("Metacell_",memberships)
          names(memberships.v5) <- names(memberships)

          mc.frag.path <- AggregateFragmentFile(input_file = fragmentFiles[["ATAC"]],
                                                tmp_path =paste0(opt$outdir,"/tmp/"),
                                                output_path = outputDirMcFragment,
                                                membership = memberships.v5)

          frag.mc <- CreateFragmentObject(mc.frag.path,
                                          cells = colnames(seurat.mc.multi))

          Fragments(seurat.mc.multi[["ATAC"]]) <- frag.mc

        }
      }
    }
  }
} else {
  seurat.mc <- readRDS(opt$inputSeuratMetacell)


  seurat.mc.multi <- SCimplify_for_Seurat_v5(seurat = seurat, seurat.mc = seurat.mc, gamma = opt$gamma)

}


if (return.seurat) {
  saveRDS(seurat.mc.multi, paste0(opt$outdir,"/seurat.multiome.mc.rds"))
} else {
  write.csv(data.frame(seurat.mc.multi$membership),
            paste0(opt$outdir,"/SuperCellMemberships.csv"))
}
