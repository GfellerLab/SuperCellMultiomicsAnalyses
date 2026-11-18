# SemiSupIntegrateLayers <- function (object, method, orig.reduction = "pca", assay = NULL, 
#           features = NULL, layers = NULL, scale.layer = "scale.data", label.col = NULL,
#           ...) 
# {
#   if (is_quosure(x = method)) {
#     method <- eval(expr = quo_get_expr(quo = method), envir = quo_get_env(quo = method))
#   }
#   if (is.character(x = method)) {
#     method <- get(x = method)
#   }
#   if (!is.function(x = method)) {
#     abort(message = "'method' must be a function for integrating layers")
#   }
#   assay <- assay %||% DefaultAssay(object = object)
#   if (inherits(x = object[[assay]], what = "SCTAssay")) {
#     layers <- "data"
#     scale.layer <- "scale.data"
#     features <- features %||% SelectSCTIntegrationFeatures(object = object, 
#                                                            assay = assay)
#   }
#   else if (inherits(x = object[[assay]], what = "StdAssay")) {
#     layers <- Layers(object = object, assay = assay, search = layers %||% 
#                        "data")
#     scale.layer <- Layers(object = object, search = scale.layer)
#     features <- features %||% VariableFeatures(object = object, 
#                                                assay = assay, nfeatures = 2000L)
#   }
#   else {
#     abort(message = "'assay' must be a v5 or SCT assay")
#   }
#   if (!is.null(scale.layer)) {
#     features <- intersect(x = features, y = Features(x = object, 
#                                                      assay = assay, layer = scale.layer))
#   }
#   if (!length(x = features)) {
#     abort(message = "None of the features provided are found in this assay")
#   }
#   if (!is.null(orig.reduction)) {
#     orig.reduction <- orig.reduction %||% DefaultDimReduc(object = object, 
#                                                           assay = assay)
#     if (!orig.reduction %in% Reductions(object = object)) {
#       abort(message = paste(sQuote(x = orig.reduction), 
#                             "is not a dimensional reduction"))
#     }
#     obj.orig <- object[[orig.reduction]]
#     if (is.null(x = DefaultAssay(object = obj.orig))) {
#       DefaultAssay(object = obj.orig) <- assay
#     }
#   }
#   value <- method(object = object[[assay]], assay = assay, 
#                   orig = obj.orig, layers = layers, scale.layer = scale.layer, 
#                   features = features, labels = object[[label.col]],..)  
#   for (i in names(x = value)) {
#     object[[i]] <- value[[i]]
#   }
#   return(object)
# }



### object is directly an assay !! 

SemiSupRPCAIntegration <- function (object = NULL, assay = NULL, layers = NULL, orig = NULL, 
                                    new.reduction = "integrated.dr", reference = NULL, features = NULL, 
                                    normalization.method = c("LogNormalize", "SCT"), dims = 1:30, 
                                    k.filter = NA, scale.layer = "scale.data", dims.to.integrate = NULL, 
                                    k.weight = 100, weight.reduction = NULL, sd.weight = 1, sample.tree = NULL, 
                                    preserve.order = FALSE, verbose = TRUE,labels.df = NULL,anchor.filtering.mtx=NULL, ...) 
{
  
  op <- options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  normalization.method <- match.arg(arg = normalization.method)
  features <- features %||% SelectIntegrationFeatures5(object = object)
  assay <- assay %||% "RNA"
  layers <- layers %||% SeuratObject::Layers(object = object, search = "data")
  ncells <- sapply(X = layers, FUN = function(x) {
    ncell <- dim(object[x])[2]
    return(ncell)
  })
  if (min(ncells) < max(dims)) {
    abort(message = "At least one layer has fewer cells than dimensions specified, please lower 'dims' accordingly.")
  }
  if (normalization.method == "SCT") {
    groups <- CreateIntegrationGroups(object, layers = layers, 
                                      scale.layer = scale.layer)
    object.sct <- CreateSeuratObject(counts = object, assay = "SCT")
    object.sct$split <- groups[, 1]
    object.list <- SplitObject(object = object.sct, split.by = "split")
    object.list <- PrepSCTIntegration(object.list = object.list, 
                                      anchor.features = features)
    object.list <- lapply(X = object.list, FUN = function(x) {
      x <- RunPCA(object = x, features = features, verbose = FALSE, 
                  npcs = max(dims))
      return(x)
    })
  }
  else {
    
    object.list <- list()
    for (i in seq_along(along.with = layers)) {
      object.list[[i]] <- suppressMessages(suppressWarnings(CreateSeuratObject(counts = NULL, 
                                                                               data = object[layers[i]][features, ])))
      VariableFeatures(object = object.list[[i]]) <- features
      object.list[[i]] <- AddMetaData(object.list[[i]],metadata = labels.df[colnames(object.list[[i]]),],col.name = colnames(labels.df))
      
      object.list[[i]] <- suppressWarnings(ScaleData(object = object.list[[i]], 
                                                     verbose = FALSE))
      object.list[[i]] <- RunPCA(object = object.list[[i]], 
                                 verbose = FALSE, npcs = max(dims))
      suppressWarnings(object.list[[i]][["RNA"]]$counts <- NULL) #Actually not needed but present in the original function
    }
  }
  
  anchors <- FindIntegrationAnchors(object.list = object.list, 
                                    anchor.features = features, scale = FALSE, reduction = "rpca", 
                                    normalization.method = normalization.method, dims = dims, 
                                    k.filter = k.filter, reference = reference, verbose = verbose, 
                                    ...)
  
  if (verbose) {
    message("initial anchors:")
    print((anchors))
  }
  
  slot(object = anchors, name = "object.list") <- lapply(X = slot(object = anchors, 
                                                                  name = "object.list"), FUN = function(x) {
                                                                    suppressWarnings(expr = x <- DietSeurat(x, features = features[1:2]))
                                                                    return(x)
                                                                  })
  
  object_merged <- IntegrateEmbeddings(anchorset = anchors,
                                       reductions = orig, new.reduction.name = paste0("unsup_",new.reduction),
                                       dims.to.integrate = dims.to.integrate, k.weight = k.weight,
                                       weight.reduction = weight.reduction, sd.weight = sd.weight,
                                       sample.tree = sample.tree, preserve.order = preserve.order,
                                       verbose = verbose)
  
  output.list <- list(object_merged[[paste0("unsup_",new.reduction)]])
  names(output.list)<- paste0("unsup_",new.reduction) 
  print(output.list)
  remove(object_merged)
  gc(verbose = F)
  
  message("unsup done")
  
  ## semi sup integration with anchor filtering
  
  if (is.null(anchor.filtering.mtx)) {
    anchor.filtering.mtx <- matrix(F,ncol = length(unique(labels.df[,1])),
                                   nrow =  length(unique(labels.df[,1])))
    colnames(anchor.filtering.mtx) <- unique(labels.df[,1])
    rownames(anchor.filtering.mtx) <- unique(labels.df[,1])
    
    diag(anchor.filtering.mtx) <- T
    anchor.filtering.mtx["unknown",] <- T
    anchor.filtering.mtx[,"unknown"] <- T
  }
  
  if (verbose) {
    message("anchor filtering matrix:")
    print(anchor.filtering.mtx)
  }
  
  filtered.anchors <- subset(anchors, group.by = colnames(labels.df)[1], ident.matrix = anchor.filtering.mtx)
  
  if (verbose) {
    message("filtered anchors:")
    print((filtered.anchors))
  }
  
  remove(anchors)
  gc(verbose = F)
  
  
  
  slot(object = filtered.anchors, name = "object.list") <- lapply(X = slot(object = filtered.anchors,
                                                                           name = "object.list"), FUN = function(x) {
                                                                             suppressWarnings(expr = x <- DietSeurat(x, features = features[1:2]))
                                                                             return(x)
                                                                           })
  
  object_merged <- IntegrateEmbeddings(anchorset = filtered.anchors,
                                       reductions = orig, new.reduction.name = paste0("semisup_",new.reduction),
                                       dims.to.integrate = dims.to.integrate, k.weight = k.weight,
                                       weight.reduction = weight.reduction, sd.weight = sd.weight,
                                       sample.tree = sample.tree, preserve.order = preserve.order,
                                       verbose = verbose)
  
  output.list[[paste0("semisup_",new.reduction)]] <- object_merged[[paste0("semisup_",new.reduction)]]
  
  ##
  print(names(output.list) )
  
  return(output.list)
}





STACASIntegration <- function (object = NULL, 
                               assay = NULL, 
                               layers = NULL, 
                               orig = NULL, 
                               new.reduction = "integrated.dr", 
                               reference = NULL, features = NULL, 
                               normalization.method = c("LogNormalize"), 
                               scale.layer = "scale.data",
                               integration.mode = c("unsup_and_semisup","semisup","unsup"),
                               dims = 1:30, scale.data = FALSE,
                               k.filter = NA, 
                               dims.to.integrate = NULL, 
                               k.weight = 100, 
                               weight.reduction = NULL, 
                               sd.weight = 1, 
                               sample.tree = NULL, 
                               preserve.order = FALSE, 
                               verbose = TRUE,
                               labels.df = NULL,
                               anchor.filtering.mtx=NULL, ...) 
{
  
  op <- options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  normalization.method <- match.arg(arg = normalization.method)
  features <- features %||% SelectIntegrationFeatures5(object = object)
  assay <- assay %||% "RNA"
  layers <- layers %||% Layers(object = object, search = "data")
  ncells <- sapply(X = layers, FUN = function(x) {
    ncell <- dim(object[x])[2]
    return(ncell)
  })
  if (min(ncells) < max(dims)) {
    abort(message = "At least one layer has fewer cells than dimensions specified, please lower 'dims' accordingly.")
  }
  labels.df[,1][is.na(labels.df[,1])] <- "unknown"
  if (verbose) {
	  print(table(labels.df[,1]))
  }

  object.list <- list()
  for (i in seq_along(along.with = layers)) {
    object.list[[i]] <- suppressMessages(suppressWarnings(CreateSeuratObject(counts = NULL, 
                                                                             data = object[layers[i]][features, ])))
    VariableFeatures(object = object.list[[i]]) <- features
    object.list[[i]] <- AddMetaData(object.list[[i]],metadata = labels.df[colnames(object.list[[i]]),],col.name = colnames(labels.df))
  }
  
  
  message("Finding anchors using STACAS...")
  
  stacas.anchors <- FindAnchors.STACAS(object.list = object.list, assay = rep(assay,length(object.list)), 
                                       anchor.features = features, reference = reference, 
                                       min.sample.size = 0, max.seed.objects = length(object.list), 
                                       genesBlockList = NULL, dims = dims, 
                                       cell.labels = colnames(labels.df)[1], 
                                       scale.data = scale.data, 
                                       verbose = verbose,...)
  
  if (verbose) {
    print("initial anchors:")
    print((stacas.anchors))
  }
  
  output.list <- list()
  if (integration.mode %in% c("unsup_and_semisup","unsup")) {
    
  slot(object = stacas.anchors, name = "object.list") <- lapply(X = slot(object = stacas.anchors, 
                                                                         name = "object.list"), FUN = function(x) {
                                                                           suppressWarnings(expr = x <- DietSeurat(x, features = features[1:2]))
                                                                           return(x)
                                                                         })
  
  unsup.sample.tree <- SampleTree.STACAS(
    anchorset = stacas.anchors,
    obj.names = names(object.list),
    semisupervised = F
  )
  
  if (verbose) {
    print("Sample tree:")
    print((unsup.sample.tree))
  }
  
  object_merged <- IntegrateEmbeddings(anchorset = stacas.anchors,
                                       reductions = orig, new.reduction.name = paste0("unsup_",new.reduction),
                                       dims.to.integrate = dims.to.integrate, k.weight = k.weight,
                                       weight.reduction = weight.reduction, sd.weight = sd.weight,
                                       sample.tree = unsup.sample.tree, preserve.order = preserve.order,
                                       verbose = verbose)
  
  output.list[[paste0("unsup_",new.reduction)]] <- object_merged[[paste0("unsup_",new.reduction)]]
  

  
  remove(object_merged)
  gc(verbose = F)
  message("unsupervised integration done")
  
  }
  
  
  
  if (integration.mode %in% c("unsup_and_semisup","semisup")) {
    
  ## semi sup integration with anchor filtering
  
  if (is.null(anchor.filtering.mtx)) {
    anchor.filtering.mtx <- matrix(F,ncol = length(unique(labels.df[,1])),
                                   nrow =  length(unique(labels.df[,1])))
    colnames(anchor.filtering.mtx) <- c(unique(labels.df[,1]))
    rownames(anchor.filtering.mtx) <- c(unique(labels.df[,1]))
    
    diag(anchor.filtering.mtx) <- T
    if ("unknown" %in% unique(labels.df[,1])) {
      anchor.filtering.mtx["unknown",] <- T
      anchor.filtering.mtx[,"unknown"] <- T
    }
  }
  
  if (verbose) {
    message("anchor filtering matrix:")
    print(anchor.filtering.mtx)
  }
  
  filtered.stacas.anchors <- subset(stacas.anchors, group.by = colnames(labels.df)[1], ident.matrix = anchor.filtering.mtx)
  
  if (verbose) {
    print("filtered anchors:")
    print((filtered.stacas.anchors))
  }
  
  remove(stacas.anchors)
  gc(verbose = F)
  
  
  
  slot(object = filtered.stacas.anchors, name = "object.list") <- lapply(X = slot(object = filtered.stacas.anchors,
                                                                                  name = "object.list"), FUN = function(x) {
                                                                                    suppressWarnings(expr = x <- DietSeurat(x, features = features[1:2]))
                                                                                    return(x)
                                                                                  })
  
  semisup.sample.tree <- SampleTree.STACAS(
    anchorset = filtered.stacas.anchors,
    obj.names = names(object.list),
    semisupervised = T
  )
  
  object_merged <- IntegrateEmbeddings(anchorset = filtered.stacas.anchors,
                                       reductions = orig, new.reduction.name = paste0("semisup_",new.reduction),
                                       dims.to.integrate = dims.to.integrate, k.weight = k.weight,
                                       weight.reduction = weight.reduction, sd.weight = sd.weight,
                                       sample.tree = semisup.sample.tree, preserve.order = preserve.order,
                                       verbose = verbose)
  
  output.list[[paste0("semisup_",new.reduction)]] <- object_merged[[paste0("semisup_",new.reduction)]]
  }
  
  return(output.list)
}

merge_atac <- function(
    sample_ids,
    features,
    assay = "ATAC",
    genome = NULL,
    annotation = NULL,
    sep = c("-", "-"),
    cache_threads = NULL,
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
