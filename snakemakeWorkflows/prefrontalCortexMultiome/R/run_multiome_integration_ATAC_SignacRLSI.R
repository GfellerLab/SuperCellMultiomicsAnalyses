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
  "input-metacell-dir",  "d", 1, "character", "Path to metacell dir with saved seurat_obj.qs2 objects to integrate [required]",
  "keep_peaks",          "a", 1, "character", "Path to peaks (.qs2) to keep to to integrate [required]",    
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
  "label-col",           "l", 1, "character", "Metadata column with labels for supervised SignacRLSI; if omitted, unsupervised STACAS is run",

  # Sample selection
  "samples",             "S", 1, "character", "Samples to process in an r vector (.qs2)",
  "cache-dir",           "C", 1, "character", "Directory where per-sample metacell objects are stored on disk [default: <output-dir>/sample_cache]",
  "cache-format",        "F", 1, "character", "Format of the cached objects: auto, qs2, qs or rds. auto uses qs2, then qs, then rds, whichever is installed [default: auto]",
  "cache-threads",       "T", 1, "integer",   "Threads used to write/read qs2/qs cache files [default: 1]",
  "reuse-cache",         "u", 0, "logical",   "Reuse metacell objects already present in the cache dir instead of recomputing them (resume a run)",
  "clean-cache",         "z", 0, "logical",   "Delete the cache dir when the run finishes successfully",
  "min-peak-fraction",   "P", 1, "double",    "Keep only ATAC peaks detected in at least this fraction of all metacells, e.g. 0.05 (0 = keep all) [default: 0]",
  "n-core-samples",      "k", 1, "integer",   "Number of largest samples used as references for both STACAS and ATAC anchors [default: 10]",

  # ATAC parameters
  "atac-comp",           "A", 1, "character", "LSI dimensions considered for ATAC, e.g. '2:30' [default: 2:30]",
  "top-features-atac",   "t", 1, "character", "min.cutoff for Signac FindTopFeatures [default: q75]",
  "depth-cor-threshold", "D", 1, "double",    "Max |correlation| between an LSI component and sequencing depth [default: 0.5]",
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
  sample_col = "SampleID",
  batch_col = "batch",
  celltype_col = "cell_type",
  n_core_samples = 10L,
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


# ---- Validate arguments -----------------------------------------------------
if (is.null(opt$input_metacell_dir)) {
  die("--input-metacell-dir is required (see --help).")
}

if (is.null(opt$keep_peaks)) {
  die("--keep-peaks is required (see --help).")
}

for (f in c(opt$functions)) {
  if (!file.exists(f)) die("File not found: ", f)
}

cache_dir <- opt$input_metacell_dir
atac_comp  <- parse_dims(opt$atac_comp)
supervised <- !is.null(opt$label_col)

sample_ids <- qs2::qs_read(opt$samples, nthreads = opt$cache_threads)

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
plot_dir <- file.path(opt$output_dir, "plots")
if (!opt$no_plots) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(opt$seed)

# ---- Libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SuperCell)
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

sample_ids  <-qs2::qs_read(opt$samples)
# Optional peak filter, decided globally from the per-sample detection counts
keep_peaks <-qs2::qs_read(opt$keep_peaks)

# ---- 5. ATAC: unintegrated Cusanovich LSI ----------------------------------
# keep_peaks is still needed in step 6 to rebuild the per-sample objects from the cache
merged.objs.atac <- time_step("Merging ATAC", merge_atac(sample_ids, keep_peaks))
log_mem("ATAC merge")
DefaultAssay(merged.objs.atac) <- "ATAC"

if (!opt$sample_col %in% colnames(merged.objs.atac[[]])) {
  die("Sample column '", opt$sample_col, "' not found in merged object metadata.")
}

                                
# Largest samples act as references for ATAC anchors
smp.core.atlas <- names(sort(table(merged.objs.atac[[opt$sample_col, drop = TRUE]]), decreasing = TRUE))
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

supervised <- !is.null(opt$label_col)
if (supervised) {
  reduc_name   <- "sup_SignacRLSI"
} else {
  reduc_name   <- "unsup_SignacRLSI"
}




time_step("ATAC unintegrated LSI", {
  merged.objs.atac <- FindTopFeatures(merged.objs.atac, min.cutoff = opt$top_features_atac)
  merged.objs.atac <- RunTFIDF(merged.objs.atac, method = 2)
  merged.objs.atac <- RunSVD(merged.objs.atac, reduction.name = "lsi",n=max(atac_comp))
})

depth.cor.cusanovich <- DepthCor(merged.objs.atac, n = max(atac_comp), reduction = "lsi")
dims.lsi.cusanovich <- intersect(
  which(abs(depth.cor.cusanovich$data$nCount_ATAC) < opt$depth_cor_threshold),
  atac_comp
)
log_msg("LSI dims retained (merged): ", paste(dims.lsi.cusanovich, collapse = ", "))
if (length(dims.lsi.cusanovich) == 0) die("No LSI dimension passed the depth-correlation filter.")

merged.objs.atac <- RunUMAP(merged.objs.atac, dims = dims.lsi.cusanovich,
                            reduction = "lsi", reduction.name = "lsi.umap")
umap_plot(merged.objs.atac, "lsi.umap", "atac_lsi_umap")



# ---- 6. ATAC: per-sample LSI + anchor-based integration --------------------
atac.var.features <- VariableFeatures(merged.objs.atac)
slim_atac <- !isTRUE(opt$no_slim_atac)

# Built one sample at a time from the cache (no SplitObject of the merged ATAC object).
# Equivalent to the split-based version: same peaks, same global variable features,
# per-sample TF-IDF and SVD.
atac.list <- lapply(sample_ids, function(sid) {
  x <- cache_load(sid)
  DefaultAssay(x) <- "ATAC"
  x <- DietSeurat(x, assays = "ATAC", features = keep_peaks)
  x <- RunTFIDF(x, method = 2)
  x <- RunSVD(x, features = intersect(atac.var.features, rownames(x)),n=max(atac_comp))
  if (slim_atac) {
    # Anchor finding only needs the TF-IDF data of the anchor features and the LSI
    # reduction. Dropping counts and all non-variable peaks keeps the per-sample
    # objects (which end up stored inside the anchor set) small.
    x <- DietSeurat(x, assays = DefaultAssay(x),
                    features = intersect(atac.var.features, rownames(x)),
                    layers = "data", dimreducs = "lsi")
  }
  x
})
remove(keep_peaks)

# Names must be the sample identifiers used in the metadata (they are matched against
# the core/reference samples later)
atac.names <- vapply(atac.list, function(x) as.character(unique(x[[opt$sample_col, drop = TRUE]]))[1], character(1))
if (!identical(unname(atac.names), as.character(sample_ids))) {
  die("Sample IDs in the cached ATAC objects ('", opt$sample_col, "' column) do not match the requested samples.")
}
names(atac.list) <- atac.names
invisible(gc())
log_mem("per-sample ATAC LSI")

depth.cor.list <- lapply(atac.list, function(x) {
  dc <- DepthCor(x, n = max(atac_comp), reduction = "lsi")
  abs(dc$data$nCount_ATAC)
})
depth.cor.df <- do.call(rbind, depth.cor.list)
dims.lsi.smp <- intersect(which(colMeans(depth.cor.df) < opt$depth_cor_threshold), atac_comp)
log_msg("LSI dims retained (per-sample): ", paste(dims.lsi.smp, collapse = ", "))
if (length(dims.lsi.smp) == 0) die("No per-sample LSI dimension passed the depth-correlation filter.")

log_msg("Nb of variable features for ATAC: ", length(VariableFeatures(merged.objs.atac)))

# Same core samples as used for STACAS
ref.idx <- core_reference_idx(atac.list)

anchors <- time_step("ATAC FindIntegrationAnchors", {
  FindIntegrationAnchors(
    object.list     = atac.list,
    k.anchor        = opt$k_anchor,
    anchor.features = VariableFeatures(merged.objs.atac),
    reduction       = "rlsi",
    reference       = ref.idx,
    dims            = dims.lsi.smp
  )
})
atac.sample.names <- names(atac.list)
remove(atac.list); invisible(gc())

if (supervised) {
	anchor_label_col <- opt$label_col
  if (!anchor_label_col %in% colnames(merged.objs.atac[[]])) {
    die("Anchor label column '", anchor_label_col, "' not found in metadata.")
  }
  labs <- unique(as.character(merged.objs.atac[[anchor_label_col, drop = TRUE]]))
  anchor.filtering.mtx <- matrix(FALSE, nrow = length(labs), ncol = length(labs),
                                 dimnames = list(labs, labs))
  diag(anchor.filtering.mtx) <- TRUE
  if ("unknown" %in% labs) {
    anchor.filtering.mtx["unknown", ] <- TRUE
    anchor.filtering.mtx[, "unknown"] <- TRUE
  }
  filtered.anchors <- subset(anchors, group.by = anchor_label_col, ident.matrix = anchor.filtering.mtx)
  log_msg("Anchor filtering: ", nrow(anchors@anchors), " -> ", nrow(filtered.anchors@anchors), " rows")
} else {
  filtered.anchors <- anchors
}

# IntegrateEmbeddings merges ALL objects stored in the anchor set (twice) but only needs
# their cell names and one assay name: the embeddings to integrate are injected as a
# temporary assay. Replace the stored objects by empty placeholders with the same
# cells and assay name so that no peak matrix is copied.
slim_anchorset_objects <- function(anchorset) {
  anchorset@object.list <- lapply(anchorset@object.list, function(x) {
    cells <- colnames(x)
    m <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                              dims = c(1L, length(cells)),
                              dimnames = list("placeholder", cells))
    suppressWarnings(CreateSeuratObject(
      counts = CreateAssayObject(counts = m, check.matrix = FALSE),
      assay = DefaultAssay(x)
    ))
  })
  anchorset
}

log_mem("anchor finding")
if (slim_atac) {
  filtered.anchors <- slim_anchorset_objects(filtered.anchors)
  remove(anchors); invisible(gc())   # `anchors` still references the full objects
  log_mem("slimming the anchor set")
}

# ---- Anchors per query sample and k.weight cap ----------------------------------
# With references, every non-reference (query) sample is integrated against the
# reference set. Seurat needs at least k.weight anchors AND k.weight distinct anchor
# cells for each query, otherwise it errors (or silently lowers k for that pair).
# Check this once, up front, rather than discovering it through repeated failures.
k_weight_start <- opt$k_weight
query.idx <- setdiff(seq_along(atac.sample.names), ref.idx)
if (length(query.idx) > 0) {
  anc <- as.data.frame(filtered.anchors@anchors)
  qstats <- do.call(rbind, lapply(query.idx, function(q) {
    a <- anc[anc$dataset2 == q & anc$dataset1 %in% ref.idx, , drop = FALSE]
    data.frame(sample = atac.sample.names[q], anchors = nrow(a),
               anchor_cells = length(unique(a$cell2)))
  }))
  remove(anc)
  log_msg("Anchors between each query sample and the reference set:")
  print(qstats, row.names = FALSE)

  if (any(qstats$anchors == 0)) {
    die("No anchors found between the reference set and: ",
        paste(qstats$sample[qstats$anchors == 0], collapse = ", "),
        ". These samples cannot be integrated; remove them (--samples) or increase --n-core-samples.")
  }
  k_cap <- min(qstats$anchors, qstats$anchor_cells)
  if (k_cap < k_weight_start) {
    log_msg("Lowering k.weight from ", k_weight_start, " to ", k_cap,
            " (smallest anchor count among query samples: ",
            qstats$sample[which.min(pmin(qstats$anchors, qstats$anchor_cells))], ")")
    k_weight_start <- k_cap
  }
}

# IntegrateEmbeddings, lowering k.weight if there are too few anchors/cells
integrate_with_retry <- function(anchorset, reduction, new_name, dims, k_weight) {
  while (k_weight > 0) {
    res <- tryCatch(
      withCallingHandlers(
        IntegrateEmbeddings(anchorset = anchorset, reductions = reduction,
                            new.reduction.name = new_name, dims.to.integrate = dims,
                            k.weight = k_weight),
        error = function(e) {
          # Log where the error happened (the plain message does not say which
          # Seurat function raised it). k.weight errors are expected and retried.
          if (!grepl("k.weight", conditionMessage(e), fixed = TRUE)) {
            calls <- head(sys.calls(), -2)
            log_msg("IntegrateEmbeddings failed: ", conditionMessage(e))
            log_msg("Call stack (outermost first):")
            for (cl in tail(calls, 15)) {
              log_msg("   ", substr(paste(deparse(cl), collapse = " "), 1, 160))
            }
          }
        }
      ),
      error = function(e) e
    )
    if (!inherits(res, "error")) {
      log_msg(new_name, " integration successful with k.weight = ", k_weight)
      return(res)
    }
    msg <- conditionMessage(res)
    if (!grepl("k.weight", msg, fixed = TRUE)) stop(res)   # unrelated error: do not retry

    m <- regmatches(msg, regexec("k.weight to less than ([0-9]+)", msg))[[1]]
    new_k <- if (length(m) == 2) as.numeric(m[2]) - 10 else k_weight - 10
    if (new_k >= k_weight) new_k <- k_weight - 10
    log_msg("Error with k.weight = ", k_weight, ": ", msg)
    log_msg("Retrying with k.weight = ", new_k)
    k_weight <- new_k
  }
  stop(new_name, " integration failed: no valid k.weight found.")
}

integrated.objs.atac <- time_step("ATAC IntegrateEmbeddings", {
  integrate_with_retry(
    anchorset  = filtered.anchors,
    reduction  = merged.objs.atac[["lsi"]],
    new_name   = reduc_name,
    dims       = 1:max(atac_comp),
    k_weight   = k_weight_start
  )
})
merged.objs.atac[[reduc_name]] <- integrated.objs.atac[[reduc_name]]
remove(integrated.objs.atac, filtered.anchors)
if (exists("anchors")) remove(anchors)
log_mem("ATAC IntegrateEmbeddings")

merged.objs.atac <- RunUMAP(merged.objs.atac, reduction = reduc_name, dims = dims.lsi.cusanovich,
                            reduction.name = paste0(reduc_name, ".umap"))


umap_plot(merged.objs.atac, paste0(reduc_name, ".umap"), paste0("atac_", reduc_name, "_umap"))


# ---- 8. Save ----------------------------------------------------------------
#out_rds <- file.path(opt$output_dir, paste0(opt$prefix, "_integrated.rds"))
#saveRDS(merged.objs, out_rds)
out_qs2 <- file.path(opt$output_dir, paste0(opt$prefix, "_ATAC_integrated.qs2"))
qs2::qs_save(merged.objs.atac, out_qs2, nthreads = opt$cache_threads)
                     
log_msg("Saved ATAC integrated object: ", out_qs2)

writeLines(capture.output(sessionInfo()),
           file.path(opt$output_dir, paste0(opt$prefix, "_sessionInfo.txt")))
if (opt$clean_cache) {
  unlink(cache_dir, recursive = TRUE)
  log_msg("Removed cache dir: ", cache_dir)
}
log_msg("Done.")
