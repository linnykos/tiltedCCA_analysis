# from https://github.com/satijalab/seurat/blob/master/R/dimensional_reduction.R
custom_umap <- function(object, assay = NULL,
                        umap.method = 'umap-learn',
                        n.components = 2L,
                        metric = 'euclidean',
                        n.epochs = 0L,
                        learning.rate = 1,
                        min.dist = 0.3,
                        spread = 1,
                        repulsion.strength = 1,
                        negative.sample.rate = 5L,
                        a = NULL,
                        b = NULL,
                        uwot.sgd = FALSE,
                        seed.use = 42L,
                        init = "spectral",
                        metric.kwds = NULL,
                        verbose = TRUE,
                        reduction.key = 'UMAP_'){
  if (umap.method != 'umap-learn') {
    warning(
      "Running UMAP on Graph objects is only supported using the umap-learn method",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  if (!reticulate::py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
  }
  if (!reticulate::py_module_available(module = 'numpy')) {
    stop("Cannot find numpy, please install through pip (e.g. pip install numpy).")
  }
  if (!reticulate::py_module_available(module = 'sklearn')) {
    stop("Cannot find sklearn, please install through pip (e.g. pip install scikit-learn).")
  }
  if (!reticulate::py_module_available(module = 'scipy')) {
    stop("Cannot find scipy, please install through pip (e.g. pip install scipy).")
  }
  stopifnot(inherits(object, "dgCMatrix"))
  
  if(init == "custom"){
    tmp <- RSpectra::svds(object, k = 2)
    init <- multiomicCCA:::.mult_mat_vec(tmp$u, tmp$d)
  } else {
    stopifnot(init == "spectral")
  }
 
  column_vec <- colnames(object)
  
  np <- reticulate::import("numpy", delay_load = TRUE)
  sp <- reticulate::import("scipy", delay_load = TRUE)
  sklearn <- reticulate::import("sklearn", delay_load = TRUE)
  umap <- reticulate::import("umap", delay_load = TRUE)
  diag(x = object) <- 0
  data <- object
  object <- sp$sparse$coo_matrix(arg1 = object)
  
  ab.params <- umap$umap_$find_ab_params(spread = spread, min_dist = min.dist)
  a <- a %||% ab.params[[1]]
  b <- b %||% ab.params[[2]]
  n.epochs <- n.epochs %||% 0L
  random.state <- sklearn$utils$check_random_state(seed = as.integer(x = seed.use))
  umap.args <- list(
    data = object,
    graph = object,
    n_components = n.components,
    initial_alpha = learning.rate,
    a = a,
    b = b,
    gamma = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    n_epochs = as.integer(x = n.epochs),
    random_state = random.state,
    init = init,
    metric = metric,
    metric_kwds = metric.kwds,
    verbose = verbose
  )
  if (numeric_version(x = umap$pkg_resources$get_distribution("umap-learn")$version) >=
      numeric_version(x = "0.5.0")) {
    umap.args <- c(umap.args, list(
      densmap = FALSE,
      densmap_kwds = NULL,
      output_dens = FALSE
    ))
  }
  embeddings <- do.call(what = umap$umap_$simplicial_set_embedding, args = umap.args)
  if (length(x = embeddings) == 2) {
    embeddings <- embeddings[[1]]
  }
  
  embeddings <- scale(x = embeddings, scale = FALSE)
  rownames(x = embeddings) <- column_vec
  colnames(x = embeddings) <- paste0("UMAP_", 1:n.components)
  
  embeddings
}