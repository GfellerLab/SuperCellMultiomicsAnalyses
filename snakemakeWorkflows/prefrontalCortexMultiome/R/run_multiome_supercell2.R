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
  "rna-h5ad",            "r", 1, "character", "Path to the RNA .h5ad file [required]",
  "atac-h5ad",           "a", 1, "character", "Path to the ATAC .h5ad file [required]",
  "functions",           "f", 1, "character", "R file defining process_multiome() [default: functions/tools.r]",
  "python",              "y", 1, "character", "Python binary for reticulate, e.g. a conda env python [default: reticulate's choice]",

  #SuperCell2
  "target_gamma",        "g", 1, "integer", "target gamma for metacell costructio",
  "label-col",           "l", 1, "character", "Metadata column with labels for supervised SuperCell2; if omitted, unsupervised STACAS is run",

  # Outputs
  "output-dir",          "o", 1, "character", "Output directory [default: results]",
  "prefix",              "p", 1, "character", "Prefix for output file names [default: multiome]",
  "no-plots",            "N", 0, "logical",   "Do not save UMAP plots",

  # Metadata columns
  "sample-col",          "s", 1, "character", "obs column with the sample identifier [default: SampleID]",
  "batch-col",           "b", 1, "character", "obs column with the batch, only used for logging [default: batch]",
  "celltype-col",        "c", 1, "character", "Metadata column used to colour UMAPs [default: cell_type]",
 
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
  "rna-comp",            "R", 1, "character", "PCA dimensions for RNA, e.g. '1:30' or '1,2,5:10' [default: 1:30]",
  "max-seed-objects",    "m", 1, "integer",   "STACAS max.seed.objects (ignored when references are set) [default: 10]",
  "stacas-min-sample-size", "M", 1, "integer", "Samples with fewer metacells than this cannot be integrated by STACAS and cause an error [default: 100]",
  "n-variable-features", "nvf", 1, "integer",   "Number of HVGs to select [default: 1000]",
  "exclude_Mito_hvg",    "Mt", 1, "logical",   "Filter mitochondrial protein genes from HVGs (after selection [default: TRUE]",
  "exclude_Ribo_hvg",    "Rb", 1, "logical",   "Filter ribosomal protein genes from HVGs (after selection [default: TRUE]",

  # ATAC parameters
  "atac-comp",           "A", 1, "character", "LSI dimensions considered for ATAC, e.g. '1:30' [default: 1:30]",
  "top-features-atac",   "t", 1, "character", "min.cutoff for Signac FindTopFeatures [default: q75]",
  "depth-cor-threshold", "d", 1, "double",    "Max |correlation| between an LSI component and sequencing depth [default: 0.5]",
  "k-anchor",            "K", 1, "integer",   "k.anchor for FindIntegrationAnchors [default: 5]",
  "no-slim-atac",        "x", 0, "logical",   "Disable memory-saving slimming of per-sample ATAC objects (keeps full peak matrices during integration)",
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
  target_gamma = 20,
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
  figure_height = 5,
  n_variable_features = 1000,
  exclude_Ribo_hvg = TRUE,
  exclude_Mito_hvg = TRUE
)

script_usage <- function(status = 0) {
  cat("Construct multiome metacells with SuperCell2.0.\n\n")
  cat(getopt::getopt(spec, usage = TRUE, command = "Rscript run_multiome_supercell2.R"))
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
if (is.null(opt$rna_h5ad) || is.null(opt$atac_h5ad)) {
  die("--rna-h5ad and --atac-h5ad are required (see --help).")
}
for (f in c(opt$rna_h5ad, opt$atac_h5ad, opt$functions)) {
  if (!file.exists(f)) die("File not found: ", f)
}

rna_comp  <- parse_dims(opt$rna_comp)
atac_comp <- parse_dims(opt$atac_comp)
supervised <- !is.null(opt$label_col)
anchor_label_col <- if (!is.null(opt$anchor_label_col)) opt$anchor_label_col else opt$label_col

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
  library(Signac)
  library(ggplot2)
})



source(opt$functions)


if (!exists("process_multiome")) die("process_multiome() not defined in ", opt$functions)



# ---- 1. Load data -----------------------------------------------------------
log_msg("Reading RNA:  ", opt$rna_h5ad)
adata_brain_rna <- read_h5ad(opt$rna_h5ad, backed = "r")
log_msg("Reading ATAC: ", opt$atac_h5ad)
adata_brain_atac <- read_h5ad(opt$atac_h5ad, backed = "r")

obs <- adata_brain_rna$obs
if (!opt$sample_col %in% colnames(obs)) die("Sample column '", opt$sample_col, "' not found in RNA obs.")

if (opt$batch_col %in% colnames(obs)) {
  smp_batch <- unique(obs[, c(opt$sample_col, opt$batch_col)])
  log_msg("Samples per batch:")
  print(table(smp_batch[[opt$batch_col]]))
}

#sample_ids <- as.character(unique(obs[[opt$sample_col]]))
sample_ids <- names(sort(table(obs[[opt$sample_col]]),decreasing = T))


if (!is.null(opt$samples)) {
  requested <- trimws(strsplit(opt$samples, ",")[[1]])
  missing <- setdiff(requested, sample_ids)
  if (length(missing) > 0) die("Unknown sample IDs: ", paste(missing, collapse = ", "))
  sample_ids <- requested
} else if (!is.null(opt$n_samples)) {
  sample_ids <- sample_ids[seq_len(min(opt$n_samples, length(sample_ids)))]
}
log_msg("Processing ", length(sample_ids), " sample(s)")

# ---- 2. Metacells per sample, streamed through a disk cache ---------------------
# Memory strategy: only ONE sample's metacell object is in memory while building.
# Each object is written to disk, and only small summaries are kept in RAM. The RNA
# and ATAC assays are then loaded and merged separately (never both at once), which
# also avoids the full copy of the combined object made by the previous workflow.
if (length(sample_ids) < 2) die("At least 2 samples are required for integration.")

cache_dir <- if (!is.null(opt$cache_dir)) opt$cache_dir else file.path(opt$output_dir, "sample_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

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

# Per-sample summary kept in memory: number of metacells, cell membership and, if
# peak filtering is requested, how many metacells detect each peak.
sample_stats <- function(sobj) {
  det <- NULL
  if (opt$min_peak_fraction > 0) {
    m <- methods::as(GetAssayData(sobj, assay = "ATAC", layer = "counts"), "CsparseMatrix")
    det <- tabulate(m@i + 1L, nbins = nrow(m))
    names(det) <- rownames(m)
  }
  list(ncells = ncol(sobj), membership = sobj@misc$membership, det = det)
}

stats <- setNames(vector("list", length(sample_ids)), sample_ids)
time_step("Metacell construction", {
  for (sid in sample_ids) {
    if (opt$reuse_cache && file.exists(cache_file(sid))) {
      log_msg("Reusing cached metacells for sample ", sid)
      sobj <- cache_load(sid)
    } else {
      log_msg("Building metacells for sample ", sid)
      sobj <- process_multiome(sid, adata_brain_rna, 
                               adata_brain_atac,
                               rna_dims = rna_comp,
                               atac_dims = atac_comp,
                               depth_cor_threshold = opt$depth_cor_threshold,
                               target_gamma = opt$target_gamma,
                               label = opt$label_col,
                               n_variable_features = opt$n_variable_features,
                               exclude_Ribo_hvg = opt$exclude_Ribo_hvg,
                               exclude_Mibo_hvg = opt$exclude_Mibo_hvg)
      cache_save(sobj, sid)
    }
    stats[[sid]] <- sample_stats(sobj)
    remove(sobj); invisible(gc())
  }
})
log_mem("metacell construction")
nb_samples <- length(sample_ids)
membership <- unlist(lapply(stats, function(x) x$membership))

# Optional peak filter, decided globally from the per-sample detection counts
keep_peaks <- NULL
if (opt$min_peak_fraction > 0) {
  total_metacells <- sum(vapply(stats, function(x) x$ncells, numeric(1)))
  peaks <- unique(unlist(lapply(stats, function(x) names(x$det))))
  det_all <- setNames(numeric(length(peaks)), peaks)
  for (x in stats) det_all[names(x$det)] <- det_all[names(x$det)] + x$det
  keep_peaks <- names(det_all)[det_all / total_metacells >= opt$min_peak_fraction]
  log_msg("Peak filter (detected in >= ", opt$min_peak_fraction * 100, "% of ", total_metacells,
          " metacells): keeping ", length(keep_peaks), " / ", length(peaks), " peaks")
  if (length(keep_peaks) == 0) die("No peak passes --min-peak-fraction ", opt$min_peak_fraction)
  remove(det_all, peaks)
}

# ---- 8. Save ----------------------------------------------------------------
#out_rds <- file.path(opt$output_dir, paste0(opt$prefix, "_integrated.rds"))
#saveRDS(merged.objs, out_rds)
out_qs2 <- file.path(opt$output_dir, paste0(opt$prefix, "_keep_peaks.qs2"))
qs2::qs_save(keep_peaks, out_qs2, nthreads = opt$cache_threads)

out_qs2 <- file.path(opt$output_dir, paste0(opt$prefix, "_sample_ids.qs2"))
qs2::qs_save(sample_ids, out_qs2, nthreads = opt$cache_threads)

out_qs2 <- file.path(opt$output_dir, paste0(opt$prefix, "_stats.qs2"))
qs2::qs_save(stats, out_qs2, nthreads = opt$cache_threads)

out_qs2 <- file.path(opt$output_dir, paste0(opt$prefix, "_membership.qs2"))
qs2::qs_save(membership, out_qs2, nthreads = opt$cache_threads)

                                
                     
writeLines(capture.output(sessionInfo()),
           file.path(opt$output_dir, paste0(opt$prefix, "_sessionInfo.txt")))
if (opt$clean_cache) {
  unlink(cache_dir, recursive = TRUE)
  log_msg("Removed cache dir: ", cache_dir)
}
log_msg("Done.")
