# Libraries ---------------------------------------------------------------
library(ArchR)
library(getopt)
library(Seurat)
library(Signac)

# Parameters --------------------------------------------------------------

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputSeurat',  'i', 1, "character", "REQUIRED: path to seurat object (.rds) containing at least 1 chromatin assay with fragments files",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'nWorkers', "w", 1, "numeric", "number of threads to use by ArchR",
  'genomeVersion',  'g', 1, "character", "REQUIRED: Genome version",
  'minTSS',     't',1, "numeric", 'Minimum TTS score required for each cell (default 0, to avoid any filtering)',
  'minFrags',     'm',1, "numeric", 'Minimum number of mapped fragments (default 0, to avoid any filtering)',
  'maxFrags',     'x',1, "numeric", 'Maximum number of mapped fragments (default Inf, to avoid any filtering)',
  'Upstream',     'u',1, "numeric", 'Number of basepairs upstream of the TSS (default 2000)',
  'Downstream',     'd',1, "numeric", 'Number of basepairs downstream of the TSS (default 100)',
  'ATACassay',     'a',1, "character", 'Name of the ChromatinAssay from which to extract the fragments files (default ATAC)',
  # 'outname',     'n',1, "character", 'Outdir path (default ./)',
  'removeArchrOutputs',     'r',0, "logical", 'Whether to remove te output files from ArchR or not',
  'runArchR',     'c',0, "logical", 'Whether to add ArchR gene activity or not',
  'pfm',     'p',1, "character", 'Path to pfm data'
), byrow=TRUE, ncol=5)

opt = getopt(spec)

if(is.null(opt$nWorkers)){
  opt$nWorkers = parallel::detectCores()-4
}

if(is.null(opt$minTSS)){
  opt$minTSS = 0
}

if(is.null(opt$minFrags)){
  opt$minFrags = 0
}

if(is.null(opt$maxFrags)){
  opt$maxFrags = Inf
}

if(is.null(opt$Upstream)){
  opt$Upstream = 2000
}

if(is.null(opt$Downstream)){
  opt$Downstream = 100
}

if(is.null(opt$ATACassay)){
  opt$ATACassay = "ATAC"
}

# options("scipen"=100, "digits"=4)

if(!opt$genomeVersion %in% c("hg38", "hg19", "mm10")){
  message("genome_version should one of the following: hg38, hg19, mm10")
}

# Load seurat object ------------------------------------------------------
seurat.obj <- readRDS(opt$inputSeurat)

fragments_file <- Fragments(seurat.obj[[opt$ATACassay]])
if(length(fragments_file[[1]]@path) == 0){
  stop("Please add a fragments object to your ChromatinAssay within seurat")
}


############################################################################
############################################################################
###                                                                      ###
###                        COMPUTE GA USING ARCHR                        ###
###                                                                      ###
############################################################################
############################################################################

if(opt$runArchR){
  # Compute GA using ArchR --------------------------------------------------

  # create the output directory
  dir.create(paste0(opt$outdir, "/ArchR_res/"), recursive = T)

  # Create arrow files ------------------------------------------------------
  # Choose genome version
  addArchRGenome(opt$genomeVersion)

  sampleNames <- paste0("sample", 1:length(fragments_file))

  start_path <- getwd()
  setwd(paste0(opt$outdir, "/ArchR_res/"))
  ArrowFiles <- createArrowFiles(
    inputFiles = sapply(fragments_file, function(x) x@path),
    sampleNames = sampleNames,
    validBarcodes = as.vector(sapply(fragments_file, function(x) x@cells)),
    minTSS = opt$minTSS,
    minFrags = opt$minFrags,
    maxFrags = opt$maxFrags,
    addTileMat = F,
    addGeneScoreMat = T, GeneScoreMatParams = list(geneUpstream = opt$Upstream, geneDownstream = opt$Downstream),
    threads = opt$nWorkers, force = T, cleanTmp = T
  )


  # Create project ----------------------------------------------------------
  print(getwd())
  print(ArrowFiles)
  setwd(start_path)
  ArchR_proj <- ArchRProject(
    ArrowFiles = paste0(opt$outdir, "/ArchR_res/",ArrowFiles),
    outputDirectory = paste0(opt$outdir, "/ArchR_res/"),
    copyArrows = F
  )


  # QC metrics --------------------------------------------------------------

  # Plotting QC metrics (most robust are TSS and the number of unique nuclear fragments)
  df <- getCellColData(ArchR_proj, select = c("log10(nFrags)", "TSSEnrichment"))
  p <- ggPoint(
    x = df[,1],
    y = df[,2],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
  ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

  plotPDF(p, name = "TSS-vs-Frags-filtered.pdf", ArchRProj = ArchR_proj, addDOC = FALSE)
  

  # Extract gene activity matrix --------------------------------------------

  GA_matrix = getMatrixFromProject(
    ArchRProj = ArchR_proj,
    useMatrix = "GeneScoreMatrix",
    binarize = FALSE
  )

  gene_names <- GA_matrix@elementMetadata$name
  ga_counts <- GA_matrix@assays@data$GeneScoreMatrix
  rownames(ga_counts) <- gene_names
  # ga_counts[1:5, 1:5]


  matching.barcodes <- unlist(lapply(1:length(sampleNames), function(i){
    seurat.barcodes <- names(fragments_file[[i]]@cells)
    names(seurat.barcodes) <- paste0(sampleNames[i], "#", fragments_file[[i]]@cells)
    return(seurat.barcodes)
  }))


  colnames(ga_counts) <- matching.barcodes[colnames(ga_counts)]
  # ga_counts[1:5, 1:5]

  seurat.obj[["GA_ArchR"]] <- CreateAssayObject(counts = ga_counts)

  if (!is.null(opt$removeArchrOutputs) && opt$removeArchrOutputs) {
    cat("Removing ArchR output files.\n")
    system(paste0("rm -r ", opt$outdir, "/ArchR_res/"))
  }

}


###########################################################################
###########################################################################
###                                                                     ###
###                       COMPUTE GA USING SIGNAC                       ###
###                                                                     ###
###########################################################################
###########################################################################


# Compute GA using Signac -------------------------------------------------
DefaultAssay(seurat.obj) <- "ATAC"

# genes <- read.table(opt$geneList)$x
gene.activities <- GeneActivity(object = seurat.obj, assay = "ATAC",
                                features = rownames(seurat.obj[["RNA"]]),
                                extend.upstream = opt$Upstream,
                                extend.downstream = opt$Downstream)
seurat.obj[["GA_Signac"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(seurat.obj) <- "GA_Signac"

seurat.obj <- NormalizeData(
  seurat.obj, assay = "GA_Signac",
  normalization.method = "LogNormalize", scale.factor = median(seurat.obj$nCount_GA_Signac) #10000 # to mimic the normalization performed in ArchR
)


############################################################################
############################################################################
###                                                                      ###
###                             CHROMVAR RUN                             ###
###                                                                      ###
############################################################################
############################################################################

# Run chromVAR ------------------------------------------------------------


if(opt$genomeVersion == "hg38"){
  bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
}

pwm <- readRDS(opt$pfm)
DefaultAssay(seurat.obj) <- "ATAC"
motif.matrix <- CreateMotifMatrix(features = granges(seurat.obj), pwm = pwm, genome = opt$genomeVersion, use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm)
seurat.obj <- SetAssayData(seurat.obj, assay = "ATAC", slot = 'motifs', new.data = motif.object)

BiocParallel::register(BiocParallel::MulticoreParam(opt$nWorkers))
seurat.obj <- RunChromVAR(
  object = seurat.obj,
  genome = bsgenome
)

saveRDS(seurat.obj, file = paste0(opt$outdir, "/seurat.multiome.activities.rds"))
