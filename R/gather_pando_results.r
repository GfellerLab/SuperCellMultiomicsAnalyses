
library(Pando)
library(tidyr)
library(dplyr)
library(purrr)
library("stringr")
library(ggplot2)
library(dplyr)
library(Seurat)
library(getopt)

source("R/functions/pando_functions.r")

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'data_path',  'i', 1, "character", "data path",
  'p_thresh', "p", 1, "numeric", "pval threshold",
  "outdir", "o", 1, "character", "output path"
), byrow=TRUE, ncol=5)

opt = getopt(spec)

n.target.genes <- vector()
n_regions <- vector()
n_tfs <- vector()
p_thresh <- opt$p_thresh

for(gamma in c(1, 10, 20, 50, 75, 100, 200)){
  if(gamma ==1){
    methods = c("sc")
  }else{
    methods = c("LSI1_norm", "LSI2_norm", "meanSC_norm", "aggregate", "LSI1_norm_svy", "LSI2_norm_svy", "meanSC_norm_svy")
  }
  for(m in methods){
    if(gamma == 1){ # load sc data only once
      grn_obj <- readRDS(paste0(opt$data_path, "/singlecells_analysis/grn_object_sc.rds"))
      # hvgs <- Pando::VariableFeatures(grn_obj)
      output_path <- paste0(opt$data_path, "/singlecells_analysis/")
      # sc.obj <- grn_obj@data
      # sc.obj <- DietSeurat(sc.obj, assays = c("RNA", "ATAC"))

    }else{
      grn_obj <- readRDS(paste0(opt$data_path, "/SuperCellMulti/g", gamma, "/grn_object_", m,".rds"))
      output_path <- paste0(opt$data_path, "/SuperCellMulti/g", gamma, "/")
    }
    grn_obj <- find_modules.GRNData(
      grn_obj,
      p_thresh = p_thresh,
      nvar_thresh = 2,
      min_genes_per_module = 1,
      rsq_thresh = 0.05
    )

    modules <- NetworkModules(grn_obj)

    metrics <- plot_module_metrics(grn_obj)

    metrics <- lapply(c(1:3), function(x) {
      metrics[[x]]$data$resolution = paste0("gamma.", gamma);
      return(metrics[[x]])
    })
    n.target.genes <- rbind(n.target.genes, cbind(metrics[[1]]$data, data.frame(method = m)))

    n_regions <- rbind(n_regions,  cbind(metrics[[2]]$data, data.frame(method = m)))

    n_tfs <- rbind(n_tfs,  cbind(metrics[[3]]$data, data.frame(method = m)))

    assign(x = paste0("modules_g", gamma, "_met_", m), modules)
    assign(x = paste0("coefs_g", gamma, "_met_", m), coef(grn_obj))
    saveRDS(modules, file = paste0(output_path, "/modules_g", gamma, "_met_", m, "_pval_", p_thresh, ".rds"))
    saveRDS(coef(grn_obj), file = paste0(output_path, "/coefs_g", gamma, "_met_", m, "_pval_", p_thresh, ".rds"))

  }
  rm(grn_obj)
  gc()
}
n.target.genes <- unique(n.target.genes)
n_regions <- unique(n_regions)
n_tfs <- unique(n_tfs)

save(n.target.genes, n_regions, n_tfs, file = paste0(opt$outdir, "pando_metrics_summary_p_thresh", p_thresh, ".RData"))
