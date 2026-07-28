#' Run PaCMAP (Pairwise Controlled Manifold Approximation)
#'
#' Runs PaCMAP, a method for dimensionality reduction for scRNA-seq data.
#' Constructs three kinds of pairs of points: neighbor pairs (pair_neighbors),
#' mid-near pairs (pair_MN), and further pairs (pair_FP) based on positional
#' relationship in the original space, and optimizes a low-dimensional
#' embedding accordingly. Described in Wang, Y., Huang, H., Rudin, C., &
#' Shaposhnik, Y. (2021). "Understanding how dimension reduction tools work:
#' an empirical approach to deciphering t-SNE, UMAP, TriMAP, and PaCMAP for
#' data visualization." Journal of Machine Learning Research, 22(201), 1-73.
#'
#' This wrapper uses the native R + Rcpp package \pkg{pacmapr}
#' (\url{https://github.com/williamsyy/pacmap-for-R}). No Python / conda
#' installation is required.
#'
#' @param object An object. This can be a Seurat object or a matrix-like object.
#'
#' @author Yiyang Sun, Haiyang Huang, Gaurav Rajesh Parikh
#' @references Wang, Y., Huang, H., Rudin, C., & Shaposhnik, Y. (2021).
#'
#' @examples
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunPaCMAP(object = pancreas_sub, features = Seurat::VariableFeatures(pancreas_sub))
#' DimPlot(pancreas_sub, reduction = "pacmap")
#'
#' @rdname RunPaCMAP
#' @export
RunPaCMAP <- function(object, ...) {
  if (inherits(object, "Seurat")) {
    RunPaCMAP.Seurat(object, ...)
  } else {
    RunPaCMAP.default(object, ...)
  }
}


#' @rdname RunPaCMAP
#' @method RunPaCMAP Seurat
#' @param reduction A character string specifying the reduction to be used as input. Default is "pca".
#' @param dims An integer vector specifying the dimensions to be used. Default is NULL.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param layer A character string specifying the layer name to be used. Default is "data".
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "pacmap".
#' @param reduction.key A character string specifying the prefix for the column names of the PaCMAP embeddings. Default is "PaCMAP_".
#'
#' @importFrom Seurat LogSeuratCommand DefaultAssay GetAssayData Embeddings
#' @export
RunPaCMAP.Seurat <- function(object, reduction = "pca", dims = NULL, features = NULL,
                             assay = NULL, layer = "data",
                             n_components = 2, n.neighbors = NULL, MN_ratio = 0.5, FP_ratio = 2,
                             distance_method = "euclidean",
                             lr = 1, num_iters = 250L, apply_pca = TRUE, init = "random",
                             reduction.name = "pacmap", reduction.key = "PaCMAP_",
                             n_threads = NULL,
                             verbose = TRUE, seed.use = 11L, ...) {
  if (is.null(dims) && is.null(features)) {
    stop("Please specify one of `dims` or `features`.")
  }
  if (!is.null(x = features)) {
    assay <- assay %||% DefaultAssay(object = object)

    data.use <- t(as.matrix(x = GetAssayData(object = object, layer = layer, assay = assay)[features, , drop = FALSE]))
    if (ncol(x = data.use) < n_components) {
      stop(
        "Please provide as many or more features than n_components: ",
        length(x = features),
        " features provided, ",
        n_components,
        " PaCMAP components requested",
        call. = FALSE
      )
    }
  } else if (!is.null(x = dims)) {
    if (!is.null(x = assay) && assay != DefaultAssay(object = object[[reduction]])) {
      warning("If both `assay` and `dims` are specified, the value of `assay` will get ignored.")
    }
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < n_components) {
      stop(
        "Please provide as many or more dims than n_components: ",
        length(x = dims),
        " dims provided, ",
        n_components,
        " PaCMAP components requested",
        call. = FALSE
      )
    }
  } else {
    stop("Please specify one of dims or features")
  }
  object[[reduction.name]] <- RunPaCMAP(
    object = data.use, assay = assay,
    n_components = n_components, n.neighbors = n.neighbors,
    MN_ratio = MN_ratio, FP_ratio = FP_ratio,
    distance_method = distance_method,
    lr = lr, num_iters = num_iters, apply_pca = apply_pca, init = init,
    reduction.key = reduction.key, n_threads = n_threads,
    verbose = verbose, seed.use = seed.use, ...
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}


#' @rdname RunPaCMAP
#' @method RunPaCMAP default
#' @importFrom Seurat CreateDimReducObject
#' @param n_components An integer specifying the number of PaCMAP components. Default is 2.
#' @param n.neighbors An integer specifying the number of neighbors considered in the k-Nearest Neighbor graph. Defaults to 10 for datasets with n <= 10000. For larger datasets the default is \code{round(10 + 15 * (log10(n) - 4))}.
#' @param MN_ratio A numeric value specifying the ratio of mid-near pairs to neighbor pairs. Default is 0.5.
#' @param FP_ratio A numeric value specifying the ratio of further pairs to neighbor pairs. Default is 2.
#' @param distance_method A character string specifying the distance metric to be used. One of "euclidean", "manhattan", "angular", "hamming". Default is "euclidean".
#' @param lr A numeric value specifying the Adam learning rate. Default is 1.
#' @param num_iters An integer or a length-3 integer vector giving the iterations for each of the three PaCMAP phases. A scalar is expanded to \code{c(100, 100, num_iters)}. Default is 250L, matching the original PaCMAP package (100 + 100 + 250 = 450 total iterations).
#' @param apply_pca A logical value indicating whether to apply PCA-to-100 preprocessing when the input has more than 100 features. Default is TRUE.
#' @param init A character string ("pca" or "random") or a numeric matrix used to initialize the low-dimensional embedding. Default is "random".
#' @param reduction.key A character string specifying the prefix for the column names of the PaCMAP embeddings. Default is "PaCMAP_".
#' @param n_threads Number of threads for the ANN and gradient steps. Default (\code{NULL}) uses \code{parallel::detectCores() - 1L}.
#' @param verbose A logical value indicating whether to print progress. Default is TRUE.
#' @param seed.use An integer specifying the random seed (passed to \code{random_state}). Default is 11.
#' @param ... Additional arguments forwarded to \code{pacmapr::pacmap}.
#' @export
RunPaCMAP.default <- function(object, assay = NULL,
                              n_components = 2, n.neighbors = NULL, MN_ratio = 0.5, FP_ratio = 2,
                              distance_method = "euclidean",
                              lr = 1, num_iters = 250L, apply_pca = TRUE, init = "random",
                              reduction.key = "PaCMAP_",
                              n_threads = NULL,
                              verbose = TRUE, seed.use = 11L, ...) {
  if (!requireNamespace("pacmapr", quietly = TRUE)) {
    stop("Package 'pacmapr' is required. Install from GitHub with:\n",
         "  remotes::install_github(\"williamsyy/pacmap-for-R\", subdir = \"pacmapr\")",
         call. = FALSE)
  }
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  X <- as.matrix(object)
  storage.mode(X) <- "double"

  res <- pacmapr::pacmap(
    X            = X,
    n_components = as.integer(n_components),
    n_neighbors  = if (is.null(n.neighbors)) NULL else as.integer(n.neighbors),
    MN_ratio     = MN_ratio,
    FP_ratio     = FP_ratio,
    distance     = distance_method,
    lr           = lr,
    num_iters    = num_iters,
    init         = init,
    apply_pca    = apply_pca,
    n_threads    = n_threads,
    random_state = if (is.null(seed.use)) NULL else as.integer(seed.use),
    verbose      = verbose,
    ...
  )

  embedding <- res$embedding
  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(x = embedding)))
  rownames(x = embedding) <- rownames(X)

  reduction <- CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  return(reduction)
}
