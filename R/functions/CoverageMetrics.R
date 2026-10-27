# ============================================================
# Metacell coverage metrics for a Seurat object
#
# Assumptions:
#   - PCA coordinates are stored in the Seurat "pca" reduction
#   - Each cell has a metacell ID in meta.data
#   - Metacell centroids are the mean PCA coordinates
#     of their member cells
#
# Metrics:
#   1. Raw D95
#   2. Density-normalized D95
#   3. Centroid Recall @ k
# ============================================================

library(Seurat)
library(RANN)


# ------------------------------------------------------------
# 1. Extract PCA coordinates
# ------------------------------------------------------------

get_pca_matrix <- function(
    seu,
    reduction = "pca",
    dims = 1:30
) {
  
  emb <- Embeddings(seu, reduction = reduction)
  
  missing_dims <- dims[dims > ncol(emb)]
  
  if (length(missing_dims) > 0) {
    stop(
      "Requested PC(s) not available: ",
      paste(missing_dims, collapse = ", ")
    )
  }
  
  emb[, dims, drop = FALSE]
}


# ------------------------------------------------------------
# 2. Calculate metacell centroids
#
# Each centroid = mean PCA coordinate of member cells
# ------------------------------------------------------------

get_metacell_centroids <- function(
    seu,
    metacell_col,
    reduction = "pca",
    dims = 1:30,
    drop_na = TRUE
) {
  
  # PCA coordinates
  pca <- get_pca_matrix(
    seu,
    reduction = reduction,
    dims = dims
  )
  
  # Metadata
  meta <- seu[[]]
  
  if (!metacell_col %in% colnames(meta)) {
    stop(
      "Column '", metacell_col,
      "' not found in Seurat metadata."
    )
  }
  
  metacell_id <- meta[[metacell_col]]
  
  # Remove cells with missing metacell assignment
  if (drop_na) {
    
    keep <- !is.na(metacell_id)
    
    pca <- pca[keep, , drop = FALSE]
    metacell_id <- metacell_id[keep]
    
  } else {
    
    if (any(is.na(metacell_id))) {
      stop(
        "Some cells have NA metacell assignments. ",
        "Set drop_na = TRUE or resolve these assignments."
      )
    }
  }
  
  metacell_id <- as.character(metacell_id)
  
  # Calculate mean PCA coordinates for each metacell
  centroids <- rowsum(
    pca,
    group = metacell_id,
    reorder = FALSE
  )
  
  # Number of cells per metacell
  sizes <- table(metacell_id)
  
  # Divide sums by number of cells
  centroids <- centroids / as.numeric(sizes)
  
  # Preserve metacell IDs
  rownames(centroids) <- names(sizes)
  
  return(centroids)
}


# ------------------------------------------------------------
# 3. Raw centroid coverage
#
# For every cell:
#
#   d_i = distance to closest metacell centroid
#
# Main metric:
#
#   D95 = 95th percentile(d_i)
# ------------------------------------------------------------

metacell_raw_coverage <- function(
    seu,
    metacell_col,
    reduction = "pca",
    dims = 1:30,
    quantiles = c(0.50, 0.90, 0.95, 0.99)
) {
  
  # PCA coordinates
  cells <- get_pca_matrix(
    seu,
    reduction = reduction,
    dims = dims
  )
  
  # Metacell centroids
  centroids <- get_metacell_centroids(
    seu,
    metacell_col = metacell_col,
    reduction = reduction,
    dims = dims
  )
  
  # Nearest centroid for each cell
  nn <- RANN::nn2(
    data = centroids,
    query = cells,
    k = 1
  )
  
  distances <- nn$nn.d[, 1]
  
  closest_centroid <- rownames(centroids)[
    nn$nn.idx[, 1]
  ]
  
  names(distances) <- rownames(cells)
  
  result <- list(
    
    D50 = unname(quantile(distances, 0.50)),
    D90 = unname(quantile(distances, 0.90)),
    D95 = unname(quantile(distances, 0.95)),
    D99 = unname(quantile(distances, 0.99)),
    
    distances = distances,
    
    closest_centroid = closest_centroid,
    
    centroids = centroids
  )
  
  return(result)
}


# ------------------------------------------------------------
# 4. Density-normalized coverage
#
# q_i =
#
#   distance(cell_i, closest centroid)
#   ----------------------------------
#   distance(cell_i, kth nearest cell)
#
# Main metric:
#
#   D95_norm = 95th percentile(q_i)
# ------------------------------------------------------------

metacell_density_normalized_coverage <- function(
    seu,
    metacell_col,
    reduction = "pca",
    dims = 1:30,
    k = 20,
    quantiles = c(0.50, 0.90, 0.95, 0.99)
) {
  
  cells <- get_pca_matrix(
    seu,
    reduction = reduction,
    dims = dims
  )
  
  centroids <- get_metacell_centroids(
    seu,
    metacell_col = metacell_col,
    reduction = reduction,
    dims = dims
  )
  
  n_cells <- nrow(cells)
  
  if (k >= n_cells) {
    stop(
      "k must be smaller than the number of cells."
    )
  }
  
  
  # ----------------------------------------------------------
  # Distance to nearest metacell centroid
  # ----------------------------------------------------------
  
  centroid_nn <- RANN::nn2(
    data = centroids,
    query = cells,
    k = 1
  )
  
  d_centroid <- centroid_nn$nn.d[, 1]
  
  closest_centroid <- rownames(centroids)[
    centroid_nn$nn.idx[, 1]
  ]
  
  
  # ----------------------------------------------------------
  # Local cellular scale
  #
  # k + 1 because each cell is its own nearest neighbour
  # when query == data
  # ----------------------------------------------------------
  
  cell_nn <- RANN::nn2(
    data = cells,
    query = cells,
    k = k + 1
  )
  
  d_local <- cell_nn$nn.d[, k + 1]
  
  
  # ----------------------------------------------------------
  # Density-normalized distance
  # ----------------------------------------------------------
  
  q <- d_centroid / d_local
  
  names(q) <- rownames(cells)
  
  
  result <- list(
    
    D50_norm = unname(
      quantile(q, 0.50)
    ),
    
    D90_norm = unname(
      quantile(q, 0.90)
    ),
    
    D95_norm = unname(
      quantile(q, 0.95)
    ),
    
    D99_norm = unname(
      quantile(q, 0.99)
    ),
    
    normalized_distances = q,
    
    centroid_distances = d_centroid,
    
    local_scale = d_local,
    
    closest_centroid = closest_centroid,
    
    k = k,
    
    centroids = centroids
  )
  
  return(result)
}


# ------------------------------------------------------------
# 5. Centroid Recall @ k
#
# For each cell:
#
#   Is nearest centroid closer than kth nearest cell?
#
# Recall@k =
#
#   fraction of cells where
#
#   d(cell, centroid) <= d(cell, kth cell)
# ------------------------------------------------------------

metacell_centroid_recall <- function(
    seu,
    metacell_col,
    reduction = "pca",
    dims = 1:30,
    k_values = c(5, 10, 20, 50)
) {
  
  cells <- get_pca_matrix(
    seu,
    reduction = reduction,
    dims = dims
  )
  
  centroids <- get_metacell_centroids(
    seu,
    metacell_col = metacell_col,
    reduction = reduction,
    dims = dims
  )
  
  n_cells <- nrow(cells)
  
  max_k <- max(k_values)
  
  if (max_k >= n_cells) {
    stop(
      "All k values must be smaller than the number of cells."
    )
  }
  
  
  # ----------------------------------------------------------
  # Nearest centroid
  # ----------------------------------------------------------
  
  centroid_nn <- RANN::nn2(
    data = centroids,
    query = cells,
    k = 1
  )
  
  d_centroid <- centroid_nn$nn.d[, 1]
  
  
  # ----------------------------------------------------------
  # kNN distances between cells
  # ----------------------------------------------------------
  
  cell_nn <- RANN::nn2(
    data = cells,
    query = cells,
    k = max_k + 1
  )
  
  # First neighbour is the cell itself
  cell_distances <- cell_nn$nn.d[, -1, drop = FALSE]
  
  
  # ----------------------------------------------------------
  # Recall for each k
  # ----------------------------------------------------------
  
  recall <- sapply(
    k_values,
    function(k) {
      
      d_k <- cell_distances[, k]
      
      mean(d_centroid <= d_k)
    }
  )
  
  names(recall) <- paste0(
    "Recall@", k_values
  )
  
  
  return(
    list(
      recall = recall,
      centroid_distances = d_centroid,
      k_values = k_values,
      centroids = centroids
    )
  )
}


# ------------------------------------------------------------
# 6. Convenience function:
#    calculate all three metrics at once
# ------------------------------------------------------------

metacell_coverage_metrics <- function(
    seu,
    metacell_col,
    reduction = "pca",
    dims = 1:30,
    k_density = 20,
    k_recall = c(5, 10, 20, 50)
) {
  
  # Raw coverage
  raw <- metacell_raw_coverage(
    seu = seu,
    metacell_col = metacell_col,
    reduction = reduction,
    dims = dims
  )
  
  
  # Density-normalized coverage
  density <- metacell_density_normalized_coverage(
    seu = seu,
    metacell_col = metacell_col,
    reduction = reduction,
    dims = dims,
    k = k_density
  )
  
  
  # Centroid recall
  recall <- metacell_centroid_recall(
    seu = seu,
    metacell_col = metacell_col,
    reduction = reduction,
    dims = dims,
    k_values = k_recall
  )
  
  
  # ----------------------------------------------------------
  # Summary table
  # ----------------------------------------------------------
  
  summary <- data.frame(
    
    metric = c(
      "Raw_D50",
      "Raw_D90",
      "Raw_D95",
      "Raw_D99",
      paste0("DensityNorm_D50_k", k_density),
      paste0("DensityNorm_D90_k", k_density),
      paste0("DensityNorm_D95_k", k_density),
      paste0("DensityNorm_D99_k", k_density)
    ),
    
    value = c(
      raw$D50,
      raw$D90,
      raw$D95,
      raw$D99,
      density$D50_norm,
      density$D90_norm,
      density$D95_norm,
      density$D99_norm
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  # Add recall metrics
  recall_table <- data.frame(
    metric = names(recall$recall),
    value = as.numeric(recall$recall),
    stringsAsFactors = FALSE
  )
  
  summary <- rbind(
    summary,
    recall_table
  )
  
  
  return(
    list(
      summary = summary,
      raw = raw,
      density_normalized = density,
      recall = recall
    )
  )
}