#!/usr/bin/env Rscript
# =============================================================================
# Multiome (RNA + ATAC) metacell integration pipeline  -- LOW-MEMORY VARIANT
#
# Same analysis as run_multiome_integration.R, with two additional memory savings:
#   * the RNA object is written to the cache while the ATAC part runs (reloaded at the
#     end), so it does not coexist with the ATAC TF-IDF/SVD temporaries;
#   * per-sample ATAC LSI objects are built one at a time from the cache instead of
#     splitting the merged ATAC object, which avoided a second full copy of the peak
#     matrices at the most memory-hungry moment.
#
#   1. Build metacells per sample (process_multiome() from functions/tools.r)
#   2. Merge samples
#   3. RNA : unintegrated PCA/UMAP + STACAS integration
#   4. ATAC: Cusanovich LSI + Signac-style RLSI integration (anchors + embeddings)
#   5. Weighted nearest neighbours (WNN) on both integrated modalities
#
# Example:
#   Rscript run_multiome_integration_lowmem.R \
#     --rna-h5ad  input/final_rna_data.h5ad \
#     --atac-h5ad input/final_atac_data.h5ad \
#     --output-dir results/ \
#     --python ~/.conda/envs/supercell_tuto_env/bin/python \
#     --n-samples 40
#
# Run with --help for all options.
# =============================================================================

if (!requireNamespace("getopt", quietly = TRUE)) {
  stop("Package 'getopt' is required: install.packages('getopt')")
}

# Columns: long name | short flag | argument (0 = none, 1 = required) | type | description
# (getopt requires a unique single-character short flag for every option)
spec <- matrix(byrow = TRUE, ncol = 5, c(
  "help",                "h", 0, "logical",   "Show this help message and exit",

  # Inputs
  "input-metacell-dir",   "d", 1,  "character", "Path to metacell dir with saved qs2 object to integrate [required]",
  "functions",           "f", 1, "character", "R file defining process_multiome() [default: functions/tools.r]",
  "python",              "y", 1, "character", "Python binary for reticulate, e.g. a conda env python [default: reticulate's choice]",
  "integration-r-fun",   "i", 1, "character", "integration function (R) path for RNA",

  # Outputs
  "output-dir",          "o", 1, "character", "Output directory [default: results]",
  "prefix",              "p", 1, "character", "Prefix for output file names [default: multiome]",
  "no-plots",            "N", 0, "logical",   "Do not save UMAP plots",

  # Metadata columns
  "sample-col",          "s", 1, "character", "obs column with the sample identifier [default: SampleID]",
  "celltype-col",        "c", 1, "character", "Metadata column used to colour UMAPs [default: cell_type]",
  "label-col",           "l", 1, "character", "Metadata column with labels for supervised STACAS; if omitted, unsupervised STACAS is run",

  # Sample selection
  "samples",             "S", 1, "character", "Samples to process in an r vector (.qs2)",
  "cache-dir",           "C", 1, "character", "Directory where per-sample metacell objects are stored on disk [default: <output-dir>/sample_cache]",
  "cache-format",        "F", 1, "character", "Format of the cached objects: auto, qs2, qs or rds. auto uses qs2, then qs, then rds, whichever is installed [default: auto]",
  "cache-threads",       "T", 1, "integer",   "Threads used to write/read qs2/qs cache files [default: 1]",
  "reuse-cache",         "u", 0, "logical",   "Reuse metacell objects already present in the cache dir instead of recomputing them (resume a run)",
  "clean-cache",         "z", 0, "logical",   "Delete the cache dir when the run finishes successfully",
  "min-peak-fraction",   "P", 1, "double",    "Keep only ATAC peaks detected in at least this fraction of all metacells, e.g. 0.05 (0 = keep all) [default: 0]",
  "n-core-samples",      "k", 1, "integer",   "Number of largest samples used as references for both STACAS and ATAC anchors [default: 10]",

  # RNA parameters
  "rna-comp",            "R", 1, "character", "PCA dimensions for RNA, e.g. '1:30' or '1,2,5:10' [default: 1:30]",
  "max-seed-objects",    "m", 1, "integer",   "STACAS max.seed.objects (ignored when references are set) [default: 10]",
  "stacas-min-sample-size", "M", 1, "integer", "Samples with fewer metacells than this cannot be integrated by STACAS and cause an error [default: 100]",
  "n-variable-features", "nvf", 1, "integer",   "Number of HVGs to select [default: 1000]",
  "exclude_Mito_hvg",    "Mt", 1, "logical",   "Filter mitochondrial protein genes from HVGs (after selection [default: TRUE]",
  "exclude_Ribo_hvg",    "Rb", 1, "logical",   "Filter ribosomal protein genes from HVGs (after selection [default: TRUE]",

  # Misc
  "seed",                "e", 1, "integer",   "Random seed [default: 1234]",
  "figure-width",        "W", 1, "double",    "Width (in) of saved figures [default: 12]",
  "figure-height",       "H", 1, "double",    "Height (in) of saved figures [default: 5]"
))

defaults <- list(
  functions = "functions/tools.r",
  output_dir = "results",
  prefix = "multiome",
  no_plots = FALSE,
  sample_col = "SampleID",
  celltype_col = "cell_type",
  n_core_samples = 10L,
  rna_comp = "1:30",
  max_seed_objects = 10L,
  stacas_min_sample_size = 100L,
  k_anchor = 5L,
  k_weight = 20L,
  cache_format = "auto",
  cache_threads = 1L,
  reuse_cache = FALSE,
  clean_cache = FALSE,
  seed = 1234L,
  figure_width = 12,
  figure_height = 5,
  n_variable_features = 1000,
  exclude_Ribo_hvg = TRUE,
  exclude_Mito_hvg = TRUE
)

script_usage <- function(status = 0) {
  cat("Integrate multiome metacells with STACAS (RNA).\n\n")
  cat(getopt::getopt(spec, usage = TRUE, command = "Rscript run_multiome_integration_lowmem.R"))
  quit(save = "no", status = status)
}

if (length(commandArgs(trailingOnly = TRUE)) == 0) script_usage(1)
opt <- getopt::getopt(spec)
if (!is.null(opt$help)) script_usage(0)

# "rna-h5ad" -> "rna_h5ad" so that options can be accessed with opt$rna_h5ad
names(opt) <- gsub("-", "_", names(opt))
for (n in names(defaults)) {
  if (is.null(opt[[n]])) opt[[n]] <- defaults[[n]]
}

print(opt)

# ---- Helpers ----------------------------------------------------------------
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

log_mem <- function(label) {
  g <- gc(verbose = FALSE)
  log_msg("R memory in use after ", label, ": ", round(sum(g[, 2]) / 1024, 2), " GB")
}

die <- function(...) {
  cat("ERROR: ", paste0(...), "\n", sep = "", file = stderr())
  quit(save = "no", status = 1)
}

# Parse "1:30", "1,2,5", "1,3:6" into an integer vector
parse_dims <- function(x) {
  parts <- strsplit(x, ",")[[1]]
  unlist(lapply(parts, function(p) {
    if (grepl(":", p)) {
      r <- as.integer(strsplit(p, ":")[[1]])
      r[1]:r[2]
    } else {
      as.integer(p)
    }
  }))
}

time_step <- function(label, expr) {
  t0 <- Sys.time()
  log_msg(label, " ...")
  out <- force(expr)
  log_msg(label, " done in ", format(round(difftime(Sys.time(), t0, units = "mins"), 2)))
  invisible(out)
}

source(opt$integration_r_fun)

# ---- Validate arguments -----------------------------------------------------
if (is.null(opt$input_metacell_dir)) {
  die("--input-metacell-dir is required (see --help).")
}
for (f in c(opt$functions)) {
  if (!file.exists(f)) die("File not found: ", f)
}

cache_dir <- opt$input_metacell_dir
rna_comp  <- parse_dims(opt$rna_comp)
supervised <- !is.null(opt$label_col)

sample_ids <- qs2::qs_read(opt$samples, nthreads = opt$cache_threads)

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
plot_dir <- file.path(opt$output_dir, "plots")
if (!opt$no_plots) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(opt$seed)

# ---- Libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  library(reticulate)
  if (!is.null(opt$python)) use_python(opt$python, required = TRUE)
  library(anndata)
  library(Seurat)
  library(SuperCell)
  library(ggplot2)
})

if (!requireNamespace("STACAS", quietly = TRUE)) {
  log_msg("STACAS not found, installing from GitHub")
  remotes::install_github("carmonalab/STACAS", upgrade = "never")
}
suppressPackageStartupMessages(library(STACAS))

source(opt$functions)
if (!exists("process_multiome")) die("process_multiome() not defined in ", opt$functions)

save_plot <- function(p, name) {
  if (opt$no_plots) return(invisible(NULL))
  f <- file.path(plot_dir, paste0(opt$prefix, "_", name, ".pdf"))
  ggsave(f, p, width = opt$figure_width, height = opt$figure_height)
  log_msg("Saved plot: ", f)
}

umap_plot <- function(obj, reduction, name) {
  group_cols <- intersect(c(opt$sample_col, opt$celltype_col), colnames(obj[[]]))
  save_plot(patchwork::wrap_plots(DimPlot(obj, reduction = reduction, group.by = group_cols[1])+ NoLegend(),
		       DimPlot(obj, reduction = reduction, group.by = group_cols[2])), name)
}

# Cache format: qs2/qs are much faster than RDS and compress well
cache_format <- tolower(opt$cache_format)
if (!cache_format %in% c("auto", "qs2", "qs", "rds")) die("--cache-format must be one of: auto, qs2, qs, rds")
if (cache_format == "auto") {
  cache_format <- if (requireNamespace("qs2", quietly = TRUE)) "qs2"
                  else if (requireNamespace("qs", quietly = TRUE)) "qs"
                  else "rds"
} else if (cache_format != "rds" && !requireNamespace(cache_format, quietly = TRUE)) {
  die("Package '", cache_format, "' is not installed (install.packages('", cache_format, "')) or use --cache-format rds")
}
log_msg("Cache format: ", cache_format)

cache_file <- function(sid) {
  file.path(cache_dir, paste0(gsub("[^A-Za-z0-9._-]", "_", sid), ".", cache_format))
}
cache_save <- function(obj, sid) {
  f <- cache_file(sid)
  switch(cache_format,
    qs2 = qs2::qs_save(obj, f, nthreads = opt$cache_threads),
    qs  = qs::qsave(obj, f, preset = "fast", nthreads = opt$cache_threads),
    rds = saveRDS(obj, f, compress = FALSE)
  )
  invisible(f)
}
cache_load <- function(sid) {
  f <- cache_file(sid)
  switch(cache_format,
    qs2 = qs2::qs_read(f, nthreads = opt$cache_threads),
    qs  = qs::qread(f, nthreads = opt$cache_threads),
    rds = readRDS(f)
  )
}

# ---- 1. Load data -----------------------------------------------------------

# Load one modality of every sample from the cache and merge
merge_modality <- function(assay, features = NULL) {
  objs <- lapply(sample_ids, function(sid) {
    x <- cache_load(sid)
    DefaultAssay(x) <- assay
    x <- DietSeurat(x, assays = assay, features = features)
    x
  })
  merged <- merge(objs[[1]], objs[-1])
  remove(objs); invisible(gc())
  DefaultAssay(merged) <- assay
  merged
}

merged.objs      <- time_step("Merging RNA",  merge_modality("RNA"))
log_mem("RNA merge")
# The ATAC object is by far the largest one. It is NOT merged here: it is built from the
# cache in step 5, after the RNA analysis, so that it never coexists with the RNA
# intermediates (normalisation, PCA, per-sample STACAS copies).
remove(stats); invisible(gc())

if (!opt$sample_col %in% colnames(merged.objs[[]])) {
  die("Sample column '", opt$sample_col, "' not found in merged object metadata.")
}

# Largest samples act as references for ATAC anchors
smp.core.atlas <- names(sort(table(merged.objs[[opt$sample_col, drop = TRUE]]), decreasing = TRUE))
smp.core.atlas <- smp.core.atlas[seq_len(min(opt$n_core_samples, length(smp.core.atlas)))]
log_msg("Core (reference) samples: ", paste(smp.core.atlas, collapse = ", "))

# Positions of the core samples in a per-sample object list (names = sample IDs).
# Used for BOTH the STACAS (RNA) and the anchor-based (ATAC) integrations so that
# the same samples act as references in each modality.
core_reference_idx <- function(obj_list, strict = TRUE) {
  idx <- match(smp.core.atlas, names(obj_list))
  if (anyNA(idx)) {
    msg <- paste0("Core samples missing from object list: ",
                  paste(smp.core.atlas[is.na(idx)], collapse = ", "))
    if (strict) die(msg) else log_msg("WARNING: ", msg)
    idx <- idx[!is.na(idx)]
  }
  idx
}

options(future.globals.maxSize= (opt$n_core_samples-2)*(1024^3))

# ---- 3. RNA: unintegrated analysis -----------------------------------------
time_step("RNA unintegrated analysis", {
  merged.objs <- NormalizeData(merged.objs)
  merged.objs <- FindVariableFeatures(merged.objs)
  if (opt$exclude_Ribo_hvg) {
      cat(paste0("Excluding Ribosomal protein genes from HVGs\n"))
      ribo <- STACAS:::get.blocklist(merged.objs)$Ribo
      ribo <- c(names(ribo),as.character(ribo))
      VariableFeatures(merged.objs) <- VariableFeatures(merged.objs)[!VariableFeatures(merged.objs) %in% ribo]
  }
  if (opt$exclude_Mito_hvg) {
      cat(paste0("Excluding Mitochondrial protein genes from HVGs\n"))
      mito <- STACAS:::get.blocklist(merged.objs)$Mito
      mito <- c(names(mito),as.character(mito))
      VariableFeatures(merged.objs) <- VariableFeatures(merged.objs)[!VariableFeatures(merged.objs) %in% mito]
  }
  cat(paste0("Kept HVGs for metacell identification: ", length(VariableFeatures(merged.objs)), "\n"))
  
  
  merged.objs <- ScaleData(merged.objs)
  merged.objs <- RunPCA(merged.objs, reduction.name = "pca")
  merged.objs <- RunUMAP(merged.objs, reduction = "pca", dims = rna_comp)
})
                                
umap_plot(merged.objs, "umap", paste0("rna_pca_umap"))


integration.genes <- VariableFeatures(merged.objs)
log_msg("Nb of variable features for RNA: ", length(integration.genes))

# ---- 4. RNA: STACAS integration --------------------------------------------
source(opt$integration_r_fun)

# Show progress bars even under Rscript: pbapply defaults to NO progress bar in

# non-interactive sessions, which makes long anchor searches look like a hang.

pbapply::pboptions(type = "timer")





invisible(gc())

if (supervised) {
  stacas_assay <- "semisup_STACAS"
  merged.objs[["input_label"]] <- merged.objs[[opt$label_col]]
  merged.objs <- IntegrateLayers(
  object = merged.objs,
  method = STACASIntegration,
  integration.mode = "semisup",
  orig.reduction = "pca",
  new.reduction = "STACAS",
  dims = rna_comp,
  verbose = T,
  reference = c(1:opt$n_core_samples),
  labels.df = merged.objs[["input_label"]]
  )
} else {
  stacas_assay <- "unsup_STACAS"
  merged.objs[["input_label"]] <- 'unknown'
  merged.objs <- IntegrateLayers(
  object = merged.objs,
  method = STACASIntegration,
  integration.mode = "unsup",
  orig.reduction = "pca",
  new.reduction = "STACAS",
  dims = rna_comp,
  verbose = T,
  reference = c(1:opt$n_core_samples),
  labels.df = merged.objs[["input_label"]]     #  "unknown" for unlabelled cells free to pair with every cell types in SemiSupRPCAIntegration
  )
}

merged.objs <- RunUMAP(merged.objs, reduction = stacas_assay,
                       reduction.name = paste0(stacas_assay, ".umap"), dims = rna_comp)
umap_plot(merged.objs, paste0(stacas_assay, ".umap"), paste0("rna_", stacas_assay, "_umap"))



# ---- 5. Save ----------------------------------------------------------------
#out_qs2 <- file.path(opt$output_dir, paste0(opt$prefix, "_RNA_integrated.rds"))
#saveRDS(merged.objs, out_rds)
out_qs2 <- file.path(opt$output_dir, paste0(opt$prefix, "_RNA_integrated.qs2"))
qs2::qs_save(merged.objs, out_qs2, nthreads = opt$cache_threads)
#out_qs2 <- file.path(opt$output_dir, paste0(opt$prefix, "_RNA_",stacas_assay,".qs2"))
#qs2::qs_save(merged.objs[[stacas_assay]], out_qs2, nthreads = opt$cache_threads)

#out_qs2 <- file.path(opt$output_dir, paste0(opt$prefix, "_RNA_pca_umap.qs2"))
#qs2::qs_save(merged.objs[["umap"]], out_qs2, nthreads = opt$cache_threads)
#out_qs2 <- file.path(opt$output_dir, paste0(opt$prefix, "_RNA_",stacas_assay,"_umap.qs2"))
#qs2::qs_save(merged.objs[[ paste0(stacas_assay, ".umap")]], out_qs2, nthreads = opt$cache_threads)
                     
log_msg("Saved integrated object: ", out_qs2)

writeLines(capture.output(sessionInfo()),
           file.path(opt$output_dir, paste0(opt$prefix, "_sessionInfo.txt")))

log_msg("Done.")
