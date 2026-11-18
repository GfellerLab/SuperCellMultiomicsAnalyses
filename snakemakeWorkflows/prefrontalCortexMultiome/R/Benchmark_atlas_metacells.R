library(Seurat)
library(Signac)
library(dplyr)
library(parallel)
library(getopt)
library(SuperCell)
library(scIntegrationMetrics) 




spec = matrix(c(
  'help',             'h', 0, "logical",   "Help about the program",
  'input',            'i', 1, "character", "REQUIRED : seurat object (.qs2)",
  "memberships",      "m", 1, "character", "REQUIRED : path to metacell memberships (.qs2)",
  "anno",             "a", 1, "character", "path to single cells h5ad file with annotation",
  'outdir',           'o', 1, "character", 'Outdir path (default ./)',
  'cell_type_col',    'c', 1, "character", 'cell_type col (default cell_type)',
  'batch_col',        'b', 1, "character", 'batch_col (default SampleID)',
  'reduction',        'r', 1, "character", "input reduction to evaluate",
  'assay',            'A', 1, "character", "assay used for saving results (ATAC will trigger depth cor check)",
  'python',           'p', 1, "character", "python path",
  "dims",             "d", 1, "character", "range of dims to consider for metacell integration evaluation (default 1:30)"
), byrow=TRUE, ncol=5)


opt = getopt(spec)

library(reticulate)
use_python(opt$python)
library(anndata)


if (is.null(opt$cell_type_col)) {
  opt$cell_type_col <- "cell_type"
}

if (is.null(opt$batch_col)) {
  opt$batch_col <- "SampleID"
}



if (!is.null(opt$dims)) {
  ci <- as.numeric(strsplit(opt$dims,split = ":")[[1]][1])
  cf <- as.numeric(strsplit(opt$dims,split = ":")[[1]][2])
  opt$dims <- c(ci:cf)
} else {
  opt$dims <- 1:30
}

check_depth_cor <- FALSE

print(opt)

merged.objs <- qs2::qs_read(opt$input)
meta.data.sc.all <- read_h5ad(opt$anno,backed = "r")$obs
cell_type_labels <- unique(meta.data.sc.all$cell_type)
memberships <-  qs2::qs_read(opt$memberships)
names(memberships) <- stringr::str_split_fixed(names(memberships) ,"\\.",2)[,2] 
merged.objs@misc$membership <- memberships

if (opt$assay == "ATAC") {
    DefaultAssay(merged.objs) <- "ATAC"
    depth.cor <- DepthCor(merged.objs,n = max(dims),reduction = opt$reduction)
    opt$dims <- which(abs(depth.cor$data$nCount_ATAC)< 0.5)
    print(opt$dims)
}

res <- SuperCell:::getIntegrationMetricsMetacells(sobj.mc = merged.objs,
                                           bio.label.col = opt$cell_type_col,
                                           batch.label.col = opt$batch_col,
                                           metricsLabels = unique(meta.data.sc.all[[opt$cell_type_col]]),
                                           dims = opt$dims,
                                           meta.data.sc =meta.data.sc.all[names(merged.objs@misc$membership),],
                                           reduction.name = opt$reduction)

res <- data.frame(res)
res$assay <- opt$assay
res$reduction <- opt$reduction

write.csv(res,paste0(opt$outdir,"/",opt$reduction,"bench_res.csv"),quote = FALSE,row.names = F)























