#' Run LocalMAP (Local Manifold Approximation Projection)
#'
#' Performs local-preserving manifold embedding (LocalMAP) for dimensionality reduction.
#' Similar to PaCMAP, but focuses on preserving local neighborhood structure.
#'
#' @param object A Seurat object or matrix-like data.
#' @param ... Additional arguments passed to default method.
#' @author Yiyang Sun
#' @export
RunLocalMAP <- function(object, ...) {
  if (inherits(object, "Seurat")) {
    RunLocalMAP.Seurat(object, ...)
  } else {
    RunLocalMAP.default(object, ...)
  }
}

#' @rdname RunLocalMAP
#' @method RunLocalMAP Seurat
#' @param object Seurat object
#' @param reduction name of input reduction (default "pca")
#' @param dims dimensions to use (integer vector)
#' @param features features to use (character vector)
#' @param assay assay name (character)
#' @param slot assay slot (default "data")
#' @param reduction.name name for storing result (default "localmap")
#' @param reduction.key prefix for embedding columns (default "LocalMAP_")
#' @param ... Passed to default
#' @importFrom Seurat LogSeuratCommand GetAssayData DefaultAssay CreateDimReducObject
#' @export
RunLocalMAP.Seurat <- function(
    object,
    reduction = "pca",
    dims = NULL,
    features = NULL,
    assay = NULL,
    slot = "data",
    reduction.name = "localmap",
    reduction.key = "LocalMAP_",
    ...
) {
  if (is.null(dims) && is.null(features)) {
    stop("Specify one of `dims` or `features`.")
  }
  if (!is.null(features)) {
    assay <- assay %||% DefaultAssay(object)
    mat <- t(as.matrix(GetAssayData(object, slot = slot, assay = assay)[features, , drop = FALSE]))
  } else {
    mat <- Embeddings(object[[reduction]])[, dims, drop = FALSE]
  }
  dr <- RunLocalMAP.default(mat, ...)
  object[[reduction.name]] <- dr
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunLocalMAP
#' @method RunLocalMAP default
#' @param object numeric matrix or data frame (cells x features)
#' @param assay assay name (optional)
#' @param n_components number of embedding dimensions (default 2)
#' @param n_neighbors integer, number of neighbors
#' @param lr learning rate (default 1)
#' @param num_iters iterations (default 250)
#' @param init initialization: "random" or "pca" (default "random")
#' @param reduction.key prefix for column names (default "LocalMAP_")
#' @param verbose logical (default TRUE)
#' @param seed.use integer random seed (default 42)
#' @param ... additional args passed to Python LocalMAP
#' @export
RunLocalMAP.default <- function(
    object, assay = NULL,
    n_components = 2L, n_neighbors = NULL,
    MN_ratio = 0.5, FP_ratio = 2.0,
    distance_method = "euclidean",
    lr = NULL,           # not used
    num_iters = NULL,    # not used
    init = "random",
    reduction.key = "LocalMAP_",
    verbose = TRUE,
    seed.use = 42L,
    ...
) {
  if (!is.null(seed.use)) set.seed(seed.use)
  if (!reticulate::py_module_available("pacmap")) {
    stop("Please install the 'pacmap' Python package (e.g. pip install pacmap)")
  }
  lm <- reticulate::import("pacmap")
  operator <- lm$LocalMAP(
    n_components = as.integer(n_components),
    n_neighbors = if (is.null(n_neighbors)) as.integer(10) else as.integer(n_neighbors),
    MN_ratio = MN_ratio,
    FP_ratio = FP_ratio,
    distance = distance_method
  )
  embedding <- operator$fit_transform(
    as.matrix(object),
    init = init
  )
  colnames(embedding) <- paste0(reduction.key, seq_len(ncol(embedding)))
  rownames(embedding) <- rownames(object)
  dr <- CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  return(dr)
}
