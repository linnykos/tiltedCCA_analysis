.initial_curve_fit <- function(cluster_vec,
                               dimred, 
                               lineage_order){
  stopifnot(all(lineage_order %in% cluster_vec),
            length(cluster_vec) == nrow(dimred),
            is.factor(cluster_vec))
  t(sapply(lineage_order, function(cluster){
    idx <- which(cluster_vec == cluster)
    Matrix::colMeans(dimred[idx,,drop = F])
  }))
}

.extract_pseudotime <- function(dimred,
                                initial_fit,
                                stretch){ # default stretch=2
  pcurve <- princurve::project_to_curve(dimred,
                                        s = initial_fit,
                                        stretch = stretch)
  pcurve$lambda
}