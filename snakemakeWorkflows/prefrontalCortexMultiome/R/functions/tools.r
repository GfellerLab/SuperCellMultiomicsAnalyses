
# ---- Python helpers (defined once) ----
load_helpers <- function() {
  py_run_string("
import numpy as np, scipy.sparse as sp, pandas as pd, gc

def mo_mask_idx(ad, col, sample):
    # obs is always in memory, so this never touches X
    return np.flatnonzero(np.asarray(ad.obs[col] == sample))

def mo_obs_names(ad, idx):
    idx = np.atleast_1d(np.asarray(idx, dtype=np.int64))
    return [str(x) for x in ad.obs_names[idx]]

def mo_var_names(ad):
    return [str(x) for x in ad.var_names]

def mo_load_T(ad, idx, order, layer):
    # read only the selected rows (idx must be sorted), return features x cells CSC float64
    idx = np.atleast_1d(np.asarray(idx, dtype=np.int64))
    sub = ad[idx]
    X = sub.X if layer is None else sub.layers[layer]
    X = sp.csr_matrix(X)
    if order is not None:
        X = X[np.atleast_1d(np.asarray(order, dtype=np.int64))]
    Xt = X.T.tocsc(copy=False)                   # zero-copy transpose
    Xt = Xt.astype('float64', copy=False)        # R sparse matrices need doubles
    Xt.indices = Xt.indices.astype(np.int32, copy=False)
    Xt.indptr  = Xt.indptr.astype(np.int32, copy=False)
    del sub, X
    gc.collect()
    return Xt

def mo_obs_df(ad, idx, cols):
    idx = np.atleast_1d(np.asarray(idx, dtype=np.int64))
    obs = ad.obs.iloc[idx]
    cols = list(obs.columns) if cols is None else list(cols)
    out = {}
    for c in cols:
        s = obs[c]
        if pd.api.types.is_bool_dtype(s):
            out[c] = s.astype(object).astype(str).tolist()
        elif pd.api.types.is_numeric_dtype(s):
            out[c] = s.astype('float64').tolist()
        else:
            out[c] = s.astype(object).where(s.notna(), '__NA__').astype(str).tolist()
    return out
")
  import_main(convert = TRUE)
}

get_py <- function(adata) {
  p <- adata$.__enclos_env__$private$.anndata
  if (is.null(p)) stop("Underlying Python AnnData not found; check ls(adata$.__enclos_env__$private)")
  p
}

obs_to_df <- function(lst, cells) {
  df <- as.data.frame(lst, stringsAsFactors = FALSE, check.names = FALSE)
  df[] <- lapply(df, function(v) if (is.character(v)) { v[v == "__NA__"] <- NA; v } else v)
  rownames(df) <- cells
  df
}

# ---- main function ----
build_multiome <- function(smp_id,
                           adata_rna, adata_atac,          # backed AnnData objects
                           rna_col  = "SampleID",
                           atac_col = "sample_id",
                           rna_layer = NULL, atac_layer = NULL,   # NULL = X
                           rna_obs_cols = NULL,                    # NULL = all columns
                           genome = NULL) {

  main    <- load_helpers()
  py_rna  <- get_py(adata_rna)
  py_atac <- get_py(adata_atac)

  # 1. indices and names from obs only (no matrix data read yet)
  idx_rna  <- main$mo_mask_idx(py_rna,  rna_col,  smp_id)
  idx_atac <- main$mo_mask_idx(py_atac, atac_col, smp_id)
  if (!length(idx_rna) || !length(idx_atac)) stop("Sample '", smp_id, "' not found in one modality")

  rna_cells  <- as.character(unlist(main$mo_obs_names(py_rna,  idx_rna)))
  atac_cells <- as.character(unlist(main$mo_obs_names(py_atac, idx_atac)))

  smp_pat    <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", smp_id)
  atac_cells <- sub(paste0("_?", smp_pat, "$"), "", atac_cells)

  if (anyDuplicated(rna_cells) || anyDuplicated(atac_cells)) stop("Duplicated cell names")

  # 2. keep only shared cells, so nothing extra is ever loaded
  common <- intersect(rna_cells, atac_cells)
  message(sprintf("RNA: %d | ATAC: %d | shared: %d", length(rna_cells), length(atac_cells), length(common)))
  if (!length(common)) stop("No shared cells; compare head(rna_cells) and head(atac_cells)")

  k <- rna_cells  %in% common; idx_rna  <- idx_rna[k];  rna_cells  <- rna_cells[k]
  k <- atac_cells %in% common; idx_atac <- idx_atac[k]; atac_cells <- atac_cells[k]

  # reorder ATAC rows (in Python) to follow the RNA cell order, only if needed
  ord <- match(rna_cells, atac_cells) - 1L
  if (identical(ord, seq_along(ord) - 1L)) ord <- NULL

  # 3. RNA -> Seurat
  counts <- main$mo_load_T(py_rna, as.integer(idx_rna), NULL, rna_layer)
  dimnames(counts) <- list(as.character(unlist(main$mo_var_names(py_rna))), rna_cells)
  meta <- obs_to_df(main$mo_obs_df(py_rna, as.integer(idx_rna), rna_obs_cols), rna_cells)

  sobj <- CreateSeuratObject(counts = counts, meta.data = meta)
  rm(counts, meta); invisible(gc())

  # 4. ATAC -> ChromatinAssay
  counts_atac <- main$mo_load_T(py_atac, as.integer(idx_atac), ord, atac_layer)
  peaks <- sub(":", "-", as.character(unlist(main$mo_var_names(py_atac))))
  dimnames(counts_atac) <- list(peaks, rna_cells)

  sobj[["ATAC"]] <- CreateChromatinAssay(counts = counts_atac, sep = c("-", "-"), genome = genome)
  rm(counts_atac); py_run_string("import gc; gc.collect()"); invisible(gc())

  sobj
}



process_multiome <- function(smp_id,
                             adata_rna, adata_atac,            # backed AnnData objects
                             # --- loading ---
                             rna_col = "SampleID",
                             atac_col = "sample_id",
                             rna_layer = NULL, atac_layer = NULL,
                             rna_obs_cols = NULL,
                             genome = NULL,
                             # --- RNA ---
                             n_variable_features = 1000,
                             exclude_Ribo_hvg = TRUE,
                             exlude_Mito_hvg = TRUE,
                             rna_dims = 1:30,
                             # --- ATAC ---
                             atac_min_cutoff = 1,
                             atac_dims = 1:30, 
                             depth_cor_threshold = 0.5, # skip LSI above
                             # --- metacells ---
                             target_gamma = 20,
                             min_n_metacells = 200,
                             label = NULL,
                             # --- output ---
                             save_path = NULL,                  # e.g. "sobj_UMARY-871.rds"
                             verbose = FALSE) {

  # 1. load only this sample's cells from both backed objects
  sobj <- build_multiome(smp_id, adata_rna, adata_atac,
                         rna_col = rna_col, atac_col = atac_col,
                         rna_layer = rna_layer, atac_layer = atac_layer,
                         rna_obs_cols = rna_obs_cols, genome = genome)

  # 2. RNA: normalize -> HVG -> scale -> PCA
  DefaultAssay(sobj) <- "RNA"
  sobj <- NormalizeData(sobj, verbose = verbose)
  sobj <- FindVariableFeatures(sobj, nfeatures = n_variable_features, verbose = verbose)
  if (exclude_Ribo_hvg) {
      cat(paste0("Excluding Ribosomal protein genes from HVGs"))
      ribo <- STACAS:::get.blocklist(merged.objs)$Ribo
      ribo <- c(names(ribo),as.character(ribo))
      VariableFeatures(sobj) <- VariableFeatures(sobj)[!VariableFeatures(sobj) %in% ribo]
  }
  if (exlude_Mito_hvg) {
      cat(paste0("Excluding Mitochondrial protein genes from HVGs"))
      mito <- STACAS:::get.blocklist(merged.objs)$Mito
      mito <- c(names(mito),as.character(mito))
      VariableFeatures(sobj) <- VariableFeatures(sobj)[!VariableFeatures(sobj) %in% mito]
  }
  cat(paste0("Kept HVGs for metacell identification: ", length(VariableFeatures(sobj))))
  sobj <- ScaleData(sobj, verbose = verbose)                    # variable features only (default)
  sobj <- RunPCA(sobj, npcs = max(rna_dims), verbose = verbose)

  # 3. ATAC: TF-IDF -> top features -> SVD (LSI)
  DefaultAssay(sobj) <- "ATAC"
  #print(sobj)
  sobj <- FindTopFeatures(sobj, min.cutoff = atac_min_cutoff, verbose = verbose)
  #print(sobj)
  sobj <- RunTFIDF(sobj, verbose = verbose)
  sobj <- RunSVD(sobj, n = max(atac_dims), verbose = verbose)

  
  depth.cor <- DepthCor(sobj, n = max(atac_dims), reduction = "lsi")
  dims.lsi <- intersect(
      which(abs(depth.cor$data$nCount_ATAC) < depth_cor_threshold),
      atac_comp
    )
  log_msg("LSI dims retained: ", paste(dims.lsi, collapse = ", "))  

  # 4. metacells
  gamma <- min(target_gamma,length(Cells(sobj))/min_n_metacells)
  sobj <- SCimplify_for_Seurat(
    sobj,
    label = label,
    assay     = c("RNA", "ATAC"),
    reduction = list("pca", "lsi"),
    dims      = list(rna_dims, dims.lsi),
    gamma     = gamma
  )

  # 5. optional save, then free memory
  if (!is.null(save_path)) saveRDS(sobj, save_path)
  invisible(gc())
  
  sobj <- RenameCells(sobj,add.cell.id = smp_id) 
  sobj@misc$new.membership <- paste0(smp_id,"_Metacell_",sobj@misc$membership)
  names(sobj@misc$new.membership) <- names(sobj@misc$membership)
  sobj@misc$membership <- sobj@misc$new.membership

  sobj
}

merge_atac <- function(
    sample_ids,
    features,
    assay = "ATAC",
    genome = NULL,
    annotation = NULL,
    sep = c("-", "-"),
    cache_threads = opt$cache_threads,
    membership = NULL
) {

  stopifnot(length(sample_ids) > 0)
  stopifnot(length(features) > 0)

  # Accept either GRanges or peak names
  if (inherits(features, "GRanges")) {
    peak_names <- names(features)
    if (is.null(peak_names)) {
      peak_names <- as.character(features)
    }
    ranges <- features
  } else {
    peak_names <- as.character(features)
    ranges <- Signac::StringToGRanges(
      peak_names,
      sep = sep
    )
  }

  # Store sparse matrices and cell metadata
  mats <- vector("list", length(sample_ids))
  metadata <- vector("list", length(sample_ids))

  for (i in seq_along(sample_ids)) {

    sid <- sample_ids[[i]]

    message("Loading ", i, "/", length(sample_ids), ": ", sid)

    x <- cache_load(sid)

    if (!assay %in% SeuratObject::Assays(x)) {
      stop("Missing assay ", assay, " in sample ", sid)
    }

    # Extract counts only
    mat <- SeuratObject::LayerData(
      x,
      assay = assay,
      layer = "counts"
    )

    # Retain requested features that exist in this sample
    keep <- intersect(peak_names, rownames(mat))
    mat <- mat[keep, , drop = FALSE]


    # Keep cell metadata only
    md <- x[[]]

    mats[[i]] <- mat
    metadata[[i]] <- md

    rm(x, mat, md)
    invisible(gc())
  }

  message("Binding sparse matrices...")

  # One final sparse matrix concatenation
  counts <- do.call(cbind, mats)

  rm(mats)
  invisible(gc())

  message("Combining cell metadata...")

  metadata <- do.call(rbind, metadata)
  metadata <- metadata[colnames(counts), , drop = FALSE]

  # Construct ChromatinAssay only once
  message("Creating ChromatinAssay...")

  chromatin_assay <- Signac::CreateChromatinAssay(
    counts = counts,
    ranges = ranges,
    genome = genome,
    annotation = annotation,
    min.cells = 0,
    min.features = 0
  )

  rm(counts)
  invisible(gc())

  # Construct final Seurat object
  message("Creating Seurat object...")

  obj <- SeuratObject::CreateSeuratObject(
    counts = chromatin_assay,
    assay = assay,
    meta.data = metadata
  )

  DefaultAssay(obj) <- assay

  if (!is.null(membership)) {
    obj@misc$membership <- membership
  }

  obj
}


