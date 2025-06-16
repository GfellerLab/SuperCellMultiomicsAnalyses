GetAssaySummary <- function (object, assay,...)
{
  UseMethod(generic = "GetAssaySummary", object = object)
}

GetAssaySummary.GRNData <- function (object, group_name, assay = NULL, verbose = TRUE)
{

  return(GetAssaySummary(object@data, group_name, assay = assay,
                         verbose = TRUE))
}

GetAssaySummary.Seurat <- function (object, group_name, assay = NULL, verbose = TRUE)
{
  print("corrected fun")
  if (is.null(assay)) {
    assay <- object@active.assay
  }
  print(group_name)
  print(assay)
  smry <- Seurat::Misc(object[[assay]])$summary[[group_name]]
  # smry <- NULL
  print(smry)
  if (is.null(smry)) {
    Pando:::log_message("Summary of \"", group_name, "\" does not yet exist.",
                        verbose = verbose)
    Pando:::log_message("Summarizing.", verbose = verbose)
    object <- aggregate_assay(object, assay = assay, group_name = group_name)
    smry <- GetAssaySummary(object, assay = assay, group_name = group_name,
                            verbose = verbose)
  }
  return(smry)
}
aggregate_assay <- function (object, group_name, fun = "mean", assay = "RNA", slot = "data")
{
  ass_mat <- Matrix::t(Seurat::GetAssayData(object, assay = assay,layer = slot))

  print(dim(ass_mat))

  groups <- as.character(object@meta.data[[group_name]])
  print(head(groups))
  agg_mat <- aggregate_matrix(ass_mat, groups = groups, fun = fun)
  print(dim(agg_mat))
  if (is.null(object@assays[[assay]]@misc$summary)) {
    object@assays[[assay]]@misc$summary <- list()
  }
  object@assays[[assay]]@misc$summary[[group_name]] <- agg_mat
  return(object)
}


infer_grn.GRNData <- function (object,
                               genes = NULL,
                               network_name = paste0(method, "_network"),
                               peak_to_gene_method = c("Signac", "GREAT"),
                               upstream = 1e+05,
                               downstream = 0,
                               extend = 1e+06,
                               only_tss = FALSE,
                               parallel = FALSE, tf_cor = 0.1, peak_cor = 0, aggregate_rna_col = NULL,
                               aggregate_peaks_col = NULL, method = c("glm", "svyglm", "glmnet", "cv.glmnet",
                                                                      "brms", "xgb", "bagging_ridge", "bayesian_ridge"), alpha = 0.5,
                               family = "gaussian", interaction_term = ":", adjust_method = "fdr",
                               scale = FALSE, verbose = TRUE, ...)
{
  method <- match.arg(method)
  peak_to_gene_method <- match.arg(peak_to_gene_method)
  object <- fit_grn_models.GRNData(object = object, genes = genes,
                           network_name = network_name, peak_to_gene_method = peak_to_gene_method,
                           upstream = upstream, downstream = downstream, extend = extend,
                           only_tss = only_tss, parallel = parallel, tf_cor = tf_cor,
                           peak_cor = peak_cor, aggregate_rna_col = aggregate_rna_col,
                           aggregate_peaks_col = aggregate_peaks_col, method = method,
                           alpha = alpha, family = family, interaction_term = interaction_term,
                           adjust_method = adjust_method, scale = scale, verbose = verbose,
                           ...)
  return(object)
}


library("sparseMatrixStats")
fit_grn_models.GRNData <- function (object, genes = NULL, network_name = paste0(method,
                                                                                "_network"), peak_to_gene_method = c("Signac", "GREAT"),
                                    upstream = 1e+05, downstream = 0, extend = 1e+06, only_tss = FALSE,
                                    peak_to_gene_domains = NULL, parallel = FALSE, tf_cor = 0.1,
                                    peak_cor = 0, aggregate_rna_col = NULL, aggregate_peaks_col = NULL,
                                    method = c("glm","svyglm", "glmnet", "cv.glmnet", "brms", "xgb", "bagging_ridge",
                                               "bayesian_ridge"), interaction_term = ":", adjust_method = "fdr",
                                    scale = FALSE, verbose = TRUE, weights = NULL, ...)
{
  method <- match.arg(method)
  peak_to_gene_method <- match.arg(peak_to_gene_method)
  Pando:::check_if_available(method)

  if(!is.null(weights)){
    aggregate_rna_col = NULL
    aggregate_peaks_col = NULL
  }
  params <- Params(object)
  motif2tf <- NetworkTFs(object)
  if (is.null(motif2tf)) {
    stop("Motif matches have not been found. Please run `find_motifs()` first.")
  }
  gene_annot <- Signac::Annotation(Pando::GetAssay(object, params$peak_assay))
  if (is.null(gene_annot)) {
    stop("Please provide a gene annotation for the ChromatinAssay.")
  }
  if (is.null(genes)) {
    genes <- VariableFeatures(object, assay = params$rna_assay)
    if (is.null(genes)) {
      stop("Please provide a set of features or run `FindVariableFeatures()`")
    }
  }
  if (is.null(aggregate_rna_col)) {
    gene_data <- Matrix::t(Pando::LayerData(object, assay = params$rna_assay, layer = "data"))
    gene_groups <- TRUE
  }else {
    if(is.numeric(object@data@meta.data[[aggregate_rna_col]])){
      object@data@meta.data[[aggregate_rna_col]] <- paste0("g",object@data@meta.data[[aggregate_rna_col]])
    }
    gene_data <- GetAssaySummary(object, assay = params$rna_assay,
                                 group_name = aggregate_rna_col, verbose = FALSE)
    gene_groups <- object@data@meta.data[[aggregate_rna_col]]
    # print(dim(gene_data))
    # print(length(gene_groups))
  }
  if (is.null(aggregate_peaks_col)) {
    peak_data <- Matrix::t(Pando::LayerData(object, assay = params$peak_assay,
                                     layer = "data"))
    peak_groups <- TRUE
  }else {
    if(is.numeric(object@data@meta.data[[aggregate_peaks_col]])){
      object@data@meta.data[[aggregate_peaks_col]] <- paste0("g",object@data@meta.data[[aggregate_peaks_col]])
    }
    peak_data <- GetAssaySummary(object, assay = params$peak_assay,
                                 group_name = aggregate_peaks_col, verbose = T)
    peak_groups <- object@data@meta.data[[aggregate_peaks_col]]
  }
  features <- intersect(gene_annot$gene_name, genes) %>% intersect(rownames(Pando::GetAssay(object, params$rna_assay)))
  gene_annot <- gene_annot[gene_annot$gene_name %in% features, ]
  regions <- NetworkRegions(object)
  peak_data <- peak_data[, regions@peaks]
  colnames(peak_data) <- rownames(regions@motifs@data)
  peaks2motif <- regions@motifs@data
  if (is.null(peak_to_gene_domains)) {
    Pando:::log_message("Selecting candidate regulatory regions near genes",
                        verbose = verbose)
    peaks_near_gene <- find_peaks_near_genes(peaks = regions@ranges,
                                             method = peak_to_gene_method, genes = gene_annot,
                                             upstream = upstream, downstream = downstream, only_tss = only_tss)
  }else {
    Pando:::log_message("Selecting candidate regulatory regions in provided domains",
                        verbose = verbose)
    peaks_near_gene <- find_peaks_near_genes(peaks = regions@ranges,
                                             method = "Signac", genes = peak_to_gene_domains,
                                             upstream = 0, downstream = 0, only_tss = FALSE)
  }

  peaks2gene <- aggregate_matrix(t(peaks_near_gene), groups = colnames(peaks_near_gene), fun = "sum")
  peaks_at_gene <- as.logical(sparseMatrixStats::colMaxs(peaks2gene))
  peaks_with_motif <- as.logical(sparseMatrixStats::rowMaxs(peaks2motif * 1))
  peaks_use <- peaks_at_gene & peaks_with_motif
  peaks2gene <- peaks2gene[, peaks_use, drop = FALSE]
  peaks2motif <- peaks2motif[peaks_use, , drop = FALSE]
  peak_data <- peak_data[, peaks_use, drop = FALSE]
  Pando:::log_message("Preparing model input", verbose = verbose)
  tfs_use <- colnames(motif2tf)
  motif2tf <- motif2tf[, tfs_use, drop = FALSE]
  Pando:::log_message("Fitting models for ", length(features), " target genes",
                      verbose = verbose)
  names(features) <- features


  model_fits <- map_par(features, function(g) {
    if (!g %in% rownames(peaks2gene)) {
      Pando:::log_message("Warning: ", g, " not found in EnsDb",
                          verbose = verbose == 2)
      return()
    }
    gene_peaks <- as.logical(peaks2gene[g, ])
    if (sum(gene_peaks) == 0) {
      Pando:::log_message("Warning: No peaks found near ", g, verbose = verbose == 2)
      return()
    }

    g_x <- gene_data[gene_groups, g, drop = FALSE]
    peak_x <- peak_data[peak_groups, gene_peaks, drop = FALSE]
    peak_x <- peak_x[,colVars(peak_x) != 0, drop = FALSE]
    if(!is.null(weights)){
      # peak_g_cor <- as(sparse_cor(peak_x, g_x, w = weights[rownames(peak_x)]), "generalMatrix")
      peak_g_cor <-  weights::wtd.cor(as.matrix(peak_x), as.matrix(g_x), weight = weights[rownames(peak_x)])
      peak_g_cor <- peak_g_cor[,1, drop = F]
      colnames(peak_g_cor) <- g
      rownames(peak_g_cor) <- colnames(peak_x)
    }else{
      peak_g_cor <- as(Pando::sparse_cor(peak_x, g_x), "generalMatrix")
    }
    peak_g_cor[is.na(peak_g_cor)] <- 0
    peaks_use <- rownames(peak_g_cor)[abs(peak_g_cor[, 1]) > peak_cor]
    if (length(peaks_use) == 0) {
      Pando:::log_message("Warning: No correlating peaks found for ", g, verbose = verbose == 2)
      return()
    }
    peak_x <- peak_x[, peaks_use, drop = FALSE]
    peak_motifs <- peaks2motif[gene_peaks, , drop = FALSE][peaks_use, , drop = FALSE]
    gene_peak_tfs <- map(rownames(peak_motifs), function(p) {
      x <- as.logical(peak_motifs[p, ])
      peak_tfs <- colMaxs(motif2tf[x, , drop = FALSE])
      peak_tfs <- colnames(motif2tf)[as.logical(peak_tfs)]
      peak_tfs <- setdiff(peak_tfs, g)
      return(peak_tfs)
    })
    names(gene_peak_tfs) <- rownames(peak_motifs)
    gene_tfs <- purrr::reduce(gene_peak_tfs, union)
    tf_x <- gene_data[gene_groups, gene_tfs, drop = FALSE]
    if(!is.null(weights)){
      # tf_g_cor <- as(sparse_cor(tf_x, g_x, w = weights[rownames(tf_x)]), "generalMatrix")
      tf_g_cor <- weights::wtd.cor(as.matrix(tf_x), as.matrix(g_x), weight = weights[rownames(tf_x)])
      tf_g_cor <- tf_g_cor[,1, drop = F]
      colnames(tf_g_cor) <- g
      rownames(tf_g_cor) <- colnames(tf_x)
    }else{
      tf_g_cor <- as(Pando::sparse_cor(tf_x, g_x), "generalMatrix")
    }
    tf_g_cor[is.na(tf_g_cor)] <- 0
    tfs_use <- rownames(tf_g_cor)[abs(tf_g_cor[, 1]) > tf_cor]
    if (length(tfs_use) == 0) {
      Pando:::log_message("Warning: No correlating TFs found for ", g, verbose = verbose == 2)
      return()
    }
    tf_g_corr_df <- as_tibble(tf_g_cor[unique(tfs_use), ,
                                       drop = F], rownames = "tf", .name_repair = "check_unique") %>%
      dplyr::rename(tf = 1, corr = 2)
    frml_string <- map(names(gene_peak_tfs), function(p) {
      peak_tfs <- gene_peak_tfs[[p]]
      peak_tfs <- peak_tfs[peak_tfs %in% tfs_use]
      if (length(peak_tfs) == 0) {
        return()
      }
      peak_name <- str_replace_all(p, "-", "_")
      tf_name <- str_replace_all(peak_tfs, "-", "_")
      formula_str <- paste(paste(peak_name, interaction_term,
                                 tf_name, sep = " "), collapse = " + ")
      return(list(tfs = peak_tfs, frml = formula_str))
    })
    frml_string <- frml_string[!map_lgl(frml_string, is.null)]
    if (length(frml_string) == 0) {
      Pando:::log_message("Warning: No valid peak:TF pairs found for ",
                          g, verbose = verbose == 2)
      return()
    }
    target <- str_replace_all(g, "-", "_")
    model_frml <- as.formula(paste0(target, " ~ ", paste0(map(frml_string,
                                                              function(x) x$frml), collapse = " + ")))
    nfeats <- sum(map_dbl(frml_string, function(x) length(x$tfs)))
    gene_tfs <- purrr::reduce(map(frml_string, function(x) x$tfs),
                              union)
    gene_x <- gene_data[gene_groups, union(g, gene_tfs), drop = FALSE]
    model_mat <- as.data.frame(cbind(gene_x, peak_x))
    if (scale)
      model_mat <- as.data.frame(scale(as.matrix(model_mat)))
    colnames(model_mat) <- str_replace_all(colnames(model_mat),
                                           "-", "_")
    Pando:::log_message("Fitting model with ", nfeats, " variables for ",
                        g, verbose = verbose == 2)
    result <- try(fit_model(model_frml, data = model_mat,
                            method = method, weights = weights[rownames(model_mat)]), silent = TRUE)
    if (any(class(result) == "try-error")) {
      Pando:::log_message("Warning: Fitting model failed for ",
                          g, verbose = verbose)
      Pando:::log_message(result, verbose = verbose == 2)
      return()
    }else {
      result$gof$nvariables <- nfeats
      result$corr <- tf_g_corr_df
      return(result)
    }
  }, verbose = verbose, parallel = parallel)
  model_fits <- model_fits[!map_lgl(model_fits, is.null)]
  if (length(model_fits) == 0) {
    Pando:::log_message("Warning: Fitting model failed for all genes.",
                        verbose = verbose)
  }
  coefs <- map_dfr(model_fits, function(x) x$coefs, .id = "target")
  coefs <- format_coefs(coefs, term = interaction_term, adjust_method = adjust_method)
  corrs <- map_dfr(model_fits, function(x) x$corr, .id = "target")
  if (nrow(coefs) > 0) {
    coefs <- suppressMessages(left_join(coefs, corrs))
  }
  gof <- map_dfr(model_fits, function(x) x$gof, .id = "target")
  params <- list()
  params[["method"]] <- method
  params[["family"]] <- family
  params[["dist"]] <- c(upstream = upstream, downstream = downstream)
  params[["only_tss"]] <- only_tss
  params[["interaction"]] <- interaction_term
  params[["tf_cor"]] <- tf_cor
  params[["peak_cor"]] <- peak_cor
  network_obj <- new(Class = "Network", features = features,
                     coefs = coefs, fit = gof, params = params)
  object@grn@networks[[network_name]] <- network_obj
  object@grn@active_network <- network_name
  return(object)
}

sparse_cor <- function (x, y = NULL, method = "pearson", allow_neg = TRUE,
                        remove_na = TRUE, remove_inf = TRUE, weights = NULL, ...)
{
  if (method == "pearson") {
    x <- Matrix(x, sparse = TRUE)
    if (!is.null(y)) {
      y <- Matrix(y, sparse = TRUE)
    }
    corr_mat <- sparse_covcor(x, y, w = weights)$cor
  }

  if (remove_na) {
    corr_mat[is.na(corr_mat)] <- 0
  }
  if (remove_inf) {
    corr_mat[is.infinite(corr_mat)] <- 1
  }
  corr_mat <- Matrix(corr_mat, sparse = TRUE)
  if (!allow_neg) {
    corr_mat[corr_mat < 0] <- 0
  }
  return(corr_mat)
}

sparse_covcor <- function (x, y = NULL, w = NULL)
{
  if (!is(x, "dgCMatrix"))
    stop("x should be a dgCMatrix")

  if(!is.null(w)){
    if (!is.numeric(w) || length(w) != nrow(x))
      stop("w should be a numeric vector of the same length as the number of rows in x")
  }

  if (!is(y, "dgCMatrix"))
    stop("y should be a dgCMatrix")
  if (nrow(x) != nrow(y))
    stop("x and y should have the same number of rows")
  # n <- nrow(x)
  muY <- colSums(y * w) / sum(w)
  # muY <- colMeans(y)
  muX <- colSums(x * w) / sum(w)
  # muX <- colMeans(x)

  wx_centered <- sweep(x, 2, muX, "-") * sqrt(w)
  wy_centered <- sweep(y, 2, muY, "-") * sqrt(w)
  covmat <- crossprod(wx_centered, wy_centered) /  (sum(w) - sum(w^2) / sum(w))
  # covmat <- (as.matrix(crossprod(x, y)) - n * tcrossprod(muX, muY))/(n - 1)
  # sdvecX <- sqrt((colSums(x^2) - n * muX^2)/(n - 1))
  # sdvecY <- sqrt((colSums(y^2) - n * muY^2)/(n - 1))

  sdvecX <- sqrt((colSums((x - muX)^2 * w)) / (sum(w) - sum(w^2) / sum(w)))
  sdvecY <- sqrt((colSums((y - muY)^2 * w)) / (sum(w) - sum(w^2) / sum(w)))
  cormat <- covmat/tcrossprod(sdvecX, sdvecY)
  return(list(cov = covmat, cor = cormat))

}

library(dplyr)
library(purrr)
library("stringr")
find_modules.GRNData <- function (object, network = DefaultNetwork(object), p_thresh = 0.05,
          rsq_thresh = 0.1, nvar_thresh = 10, min_genes_per_module = 5)
{
  params <- Params(object)
  regions <- NetworkRegions(object)
  net_obj <- GetNetwork(object, network = network)
  net_obj <- find_modules.Network(net_obj, p_thresh = p_thresh, rsq_thresh = rsq_thresh,
                          nvar_thresh = nvar_thresh, min_genes_per_module = min_genes_per_module)
  modules <- NetworkModules(net_obj)
  reg2peaks <- rownames(Pando::GetAssay(object, assay = params$peak_assay))[regions@peaks]
  names(reg2peaks) <- Signac::GRangesToString(regions@ranges)
  peaks_pos <- modules@features$regions_pos %>% map(function(x) unique(reg2peaks[x]))
  peaks_neg <- modules@features$regions_neg %>% map(function(x) unique(reg2peaks[x]))
  modules@features[["peaks_pos"]] <- peaks_pos
  modules@features[["peaks_neg"]] <- peaks_neg
  object@grn@networks[[network]]@modules <- modules
  return(object)
}

find_modules.Network <- function (object, p_thresh = 0.05, rsq_thresh = 0.1, nvar_thresh = 10,
          min_genes_per_module = 5, xgb_method = c("tf", "target"),
          xgb_top = 50, verbose = TRUE)
{
  fit_method <- NetworkParams(object)$method
  xgb_method <- match.arg(xgb_method)
  if (!fit_method %in% c("glm","svyglm", "cv.glmnet", "glmnet", "brms",
                         "xgb")) {
    stop(paste0("find_modules() is not yet implemented for \"",
                fit_method, "\" models"))
  }
  models_use <- Pando::gof(object) %>% filter(rsq > rsq_thresh &
                                         nvariables > nvar_thresh) %>% pull(target) %>% unique()
  modules <- coef(object) %>% filter(target %in% models_use)
  if (fit_method %in% c("cv.glmnet", "glmnet")) {
    modules <- modules %>% filter(estimate != 0)
  }
  else if (fit_method == "xgb") {
    modules <- modules %>% group_by_at(xgb_method) %>% top_n(xgb_top,
                                                             gain) %>% mutate(estimate = sign(corr) * gain)
  }
  else {
    modules <- modules %>% filter(ifelse(is.na(padj), T,
                                         padj < p_thresh))
  }
  modules <- modules %>% group_by(target) %>% mutate(nvars = n()) %>%
    group_by(target, tf) %>% mutate(tf_sites_per_gene = n()) %>%
    group_by(target) %>% mutate(tf_per_gene = length(unique(tf)),
                                peak_per_gene = length(unique(region))) %>% group_by(tf) %>%
    mutate(gene_per_tf = length(unique(target))) %>% group_by(target,
                                                              tf)
  if (fit_method %in% c("cv.glmnet", "glmnet","svyglm", "xgb")) {
    modules <- modules %>% reframe(estimate = sum(estimate),
                                   n_regions = peak_per_gene, n_genes = gene_per_tf,
                                   n_tfs = tf_per_gene, regions = paste(region, collapse = ";"))
  }
  else {
    modules <- modules %>% reframe(estimate = sum(estimate),
                                   n_regions = peak_per_gene, n_genes = gene_per_tf,
                                   n_tfs = tf_per_gene, regions = paste(region, collapse = ";"),
                                   pval = min(pval), padj = min(padj))
  }
  modules <- modules %>% distinct() %>% arrange(tf)
  module_pos <- modules %>% filter(estimate > 0) %>% group_by(tf) %>%
    filter(n() > min_genes_per_module) %>% group_split() %>%
    {
      names(.) <- map_chr(., function(x) x$tf[[1]])
      .
    } %>% map(function(x) x$target)
  module_neg <- modules %>% filter(estimate < 0) %>% group_by(tf) %>%
    filter(n() > min_genes_per_module) %>% group_split() %>%
    {
      names(.) <- map_chr(., function(x) x$tf[[1]])
      .
    } %>% map(function(x) x$target)
  regions_pos <- modules %>% filter(estimate > 0) %>% group_by(tf) %>%
    filter(n() > min_genes_per_module) %>% group_split() %>%
    {
      names(.) <- map_chr(., function(x) x$tf[[1]])
      .
    } %>% map(function(x) unlist(str_split(x$regions, ";")))
  regions_neg <- modules %>% filter(estimate < 0) %>% group_by(tf) %>%
    filter(n() > min_genes_per_module) %>% group_split() %>%
    {
      names(.) <- map_chr(., function(x) x$tf[[1]])
      .
    } %>% map(function(x) unlist(str_split(x$regions, ";")))
  module_feats <- list(genes_pos = module_pos, genes_neg = module_neg,
                       regions_pos = regions_pos, regions_neg = regions_neg)
  Pando:::log_message(paste0("Found ", length(unique(modules$tf)),
                     " TF modules"), verbose = verbose)
  module_meta <- select(modules, tf, target, everything())
  object@modules@meta <- module_meta
  object@modules@features <- module_feats
  object@modules@params <- list(p_thresh = p_thresh, rsq_thresh = rsq_thresh,
                                nvar_thresh = nvar_thresh, min_genes_per_module = min_genes_per_module)
  return(object)
}

# Fitting functions -------------------------------------------------------


fit_model <- function(
    formula,
    data,
    method = c('glm', 'svyglm','glmnet', 'cv.glmnet', 'brms', 'xgb', 'bagging_ridge', 'bayesian_ridge'),
    family = gaussian,
    alpha = 1,
    weights = NULL,
    ...
){
  if(is.null(weights)){
    weights <- rep(1, nrow(data))
  }
  # Match args
  method <- match.arg(method)
  result <- switch(
    method,
    'glm' = fit_glm(formula, data, family=family, weights = weights, ...),
    'svyglm' = fit_svyglm(formula, data, family=family, weights = weights, ...),
    'glmnet' = fit_glmnet(formula, data, family=family, alpha=alpha, weights = weights, ...),
    'cv.glmnet' = fit_cvglmnet(formula, data, family=family, alpha=alpha, weights = weights, ...),
    'brms' = fit_brms(formula, data, family=family, ...),
    'xgb' = fit_xgb(formula, data, ...),
    'bagging_ridge' = fit_bagging_ridge(formula, data, alpha=alpha, ...),
    'bayesian_ridge' = fit_bayesian_ridge(formula, data, ...)
  )
  return(result)
}

fit_glm <- function(formula, data, family=gaussian, weights = NULL, ...){
  fit <- suppressWarnings(glm(formula, data=data, family=family, weights = weights, ...))
  s <- summary(fit)
  gof <- tibble(
    rsq = with(s, 1 - deviance/null.deviance)
  )
  coefs <- as_tibble(s$coefficients, rownames='term')
  colnames(coefs) <- c('term', 'estimate', 'std_err', 'statistic', 'pval')
  return(list(gof=gof, coefs=coefs))
}

fit_svyglm <- function(formula, data, family=gaussian, weights = NULL, ...){
  design_obj <- survey::svydesign(ids = ~1, data = data, weights = weights)
  fit <- suppressWarnings(survey::svyglm(formula, design = design_obj, family=family, ...))
  s <- summary(fit)
  gof <- tibble(
    rsq = with(s, 1 - deviance/null.deviance)
  )
  coefs <- as_tibble(s$coefficients, rownames='term')
  colnames(coefs) <- c('term', 'estimate', 'std_err', 'statistic', 'pval')
  return(list(gof=gof, coefs=coefs))
}


fit_glmnet <- function(
    formula,
    data,
    family = gaussian,
    alpha = 0.5, weights = NULL,
    ...
){
  fit <- glmnetUtils::glmnet(
    formula,
    data = data,
    family = family,
    alpha = alpha,
    weights = weights,
    ...
  )
  class(fit) <- 'glmnet'
  which_max <- which(fit$dev.ratio > max(fit$dev.ratio) * 0.95)[1]
  lambda_choose <- fit$lambda[which_max]
  gof <- tibble(
    lambda = lambda_choose,
    rsq = fit$dev.ratio[which_max],
    alpha = alpha
  )
  coefs <- as_tibble(as.matrix(coef(fit, s=lambda_choose)), rownames='term')
  colnames(coefs) <- c('term', 'estimate')
  return(list(gof=gof, coefs=coefs))
}
