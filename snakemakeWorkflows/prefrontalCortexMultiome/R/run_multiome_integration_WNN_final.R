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
  "rna-integrated",            "r", 1, "character", "Path to the RNA .h5ad file [required]",
  "atac-integrated",           "a", 1, "character", "Path to the ATAC .h5ad file [required]",
  "functions",           "f", 1, "character", "R file defining process_multiome() [default: functions/tools.r]",
  "integration-r-fun",   "i", 1, "character", "integration function (R) path for RNA",

  # Outputs
  "output-dir",          "o", 1, "character", "Output directory [default: results]",
  "prefix",              "p", 1, "character", "Prefix for output file names [default: multiome]",
  "no-plots",            "N", 0, "logical",   "Do not save UMAP plots",

  # Metadata columns
  "sample-col",          "s", 1, "character", "obs column with the sample identifier [default: SampleID]",
  "batch-col",           "b", 1, "character", "obs column with the batch, only used for logging [default: batch]",
  "celltype-col",        "c", 1, "character", "Metadata column used to colour UMAPs [default: cell_type]",
  "label-col",           "l", 1, "character", "Metadata column with labels for supervised STACAS; if omitted, unsupervised STACAS is run",
  "anchor-label-col",    "L", 1, "character", "Metadata column with labels used to filter ATAC anchors, supervised mode only [default: --label-col]",

  # Sample selection
  "n-samples",           "n", 1, "integer",   "Process only the first N samples [default: all]",
  "samples",             "S", 1, "character", "Comma-separated sample IDs to process (overrides --n-samples)",
  "cache-dir",           "C", 1, "character", "Directory where per-sample metacell objects are stored on disk [default: <output-dir>/sample_cache]",
  "cache-format",        "F", 1, "character", "Format of the cached objects: auto, qs2, qs or rds. auto uses qs2, then qs, then rds, whichever is installed [default: auto]",
  "cache-threads",       "T", 1, "integer",   "Threads used to write/read qs2/qs cache files [default: 1]",
  "keep-rna-in-memory",  "O", 0, "logical",   "Do not offload the RNA object to the cache during the ATAC phase (faster, but uses more RAM)",
  "reuse-cache",         "u", 0, "logical",   "Reuse metacell objects already present in the cache dir instead of recomputing them (resume a run)",
  "clean-cache",         "z", 0, "logical",   "Delete the cache dir when the run finishes successfully",
  "min-peak-fraction",   "P", 1, "double",    "Keep only ATAC peaks detected in at least this fraction of all metacells, e.g. 0.05 (0 = keep all) [default: 0]",
  "n-core-samples",      "k", 1, "integer",   "Number of largest samples used as references for both STACAS and ATAC anchors [default: 10]",

  # RNA parameters
  "rna-reduc",           "x", 1, "character", "RNA reduction to use",
  "rna-comp",            "R", 1, "character", "PCA dimensions for RNA, e.g. '1:30' or '1,2,5:10' [default: 1:30]",
  "max-seed-objects",    "m", 1, "integer",   "STACAS max.seed.objects (ignored when references are set) [default: 10]",
  "stacas-min-sample-size", "M", 1, "integer", "Samples with fewer metacells than this cannot be integrated by STACAS and cause an error [default: 100]",

  # ATAC parameters
  "atac-reduc",          "y", 1, "character", "ATAC reduction to use",
  "atac-comp",           "A", 1, "character", "LSI dimensions considered for ATAC, e.g. '2:30' [default: 2:30]",
  "top-features-atac",   "t", 1, "character", "min.cutoff for Signac FindTopFeatures [default: q75]",
  "depth-cor-threshold", "d", 1, "double",    "Max |correlation| between an LSI component and sequencing depth [default: 0.5]",
  "k-anchor",            "K", 1, "integer",   "k.anchor for FindIntegrationAnchors [default: 5]",
  "k-weight",            "w", 1, "integer",   "Initial k.weight for IntegrateEmbeddings, reduced automatically on error [default: 20]",

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
  batch_col = "batch",
  celltype_col = "cell_type",
  n_core_samples = 10L,
  rna_comp = "1:30",
  max_seed_objects = 10L,
  stacas_min_sample_size = 100L,
  atac_comp = "1:30",
  top_features_atac = "q75",
  depth_cor_threshold = 0.5,
  k_anchor = 5L,
  k_weight = 20L,
  no_slim_atac = FALSE,
  cache_format = "auto",
  cache_threads = 1L,
  keep_rna_in_memory = FALSE,
  reuse_cache = FALSE,
  clean_cache = FALSE,
  min_peak_fraction = 0,
  seed = 1234L,
  figure_width = 12,
  figure_height = 5
)

script_usage <- function(status = 0) {
  cat("Integrate multiome metacells with WNN.\n\n")
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


# ---- Validate arguments -----------------------------------------------------
if (is.null(opt$rna_integrated) || is.null(opt$atac_integrated)) {
  die("--rnarna_integratedh5ad and --atacatac_integratedh5ad are required (see --help).")
}
for (f in c(opt$rna_h5ad, opt$atac_h5ad, opt$functions)) {
  if (!file.exists(f)) die("File not found: ", f)
}

rna_comp  <- parse_dims(opt$rna_comp)
atac_comp <- parse_dims(opt$atac_comp)

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
plot_dir <- file.path(opt$output_dir, "plots")
if (!opt$no_plots) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(opt$seed)

# ---- Libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(SuperCell)
  library(Signac)
  library(ggplot2)
})



source(opt$functions)

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

# ---- 1. Load data -----------------------------------------------------------
log_msg("Reading RNA:  ", opt$rna_integrated)
rna_integrated <- qs2::qs_read(opt$rna_integrated)
log_msg("Reading ATAC: ", opt$atac_integrated)
merged.objs <- qs2::qs_read(opt$atac_integrated)


# Bring RNA assay and all RNA reductions back into the main object

merged.objs[["RNA"]] <- rna_integrated[["RNA"]]
for (r in Reductions(rna_integrated)) {
  merged.objs[[r]] <- rna_integrated[[r]]
}
remove(rna_integrated); invisible(gc())

depth.cor <- DepthCor(merged.objs, n = max(atac_comp), reduction = opt$atac_reduc)
  dims.lsi <- intersect(
      which(abs(depth.cor$data$nCount_ATAC) < opt$depth_cor_threshold),
      atac_comp
    )
  log_msg("Corrected ATAC dims retained: ", paste(dims.lsi, collapse = ", "))

# ---- 7. WNN -----------------------------------------------------------------
time_step("WNN", {
  merged.objs <- FindMultiModalNeighbors(
    merged.objs,
    reduction.list = list( opt$rna_reduc,  opt$atac_reduc),
    dims.list = list(rna_comp, dims.lsi),
    modality.weight.name = "RNA.weight"
  )
  merged.objs <- RunUMAP(merged.objs, nn.name = "weighted.nn",
                         reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
})
umap_plot(merged.objs, "wnn.umap", "wnn_umap")

# ---- 8. Save ----------------------------------------------------------------
#out_rds <- file.path(opt$output_dir, paste0(opt$prefix, "_integrated.rds"))
#saveRDS(merged.objs, out_rds)
out_qs2 <- file.path(opt$output_dir, paste0(opt$prefix, "_integrated.qs2"))
qs2::qs_save(merged.objs, out_qs2, nthreads = opt$cache_threads)
                     
log_msg("Saved integrated object: ", out_qs2)

writeLines(capture.output(sessionInfo()),
           file.path(opt$output_dir, paste0(opt$prefix, "_sessionInfo.txt")))
if (opt$clean_cache) {
  unlink(cache_dir, recursive = TRUE)
  log_msg("Removed cache dir: ", cache_dir)
}
log_msg("Done.")