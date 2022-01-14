form_kernel_matrix <- function(mat, K, 
                               metacell_clustering, 
                               clustering_hierarchy,
                               symmetrize_func = c("max", "min", "average"),
                               verbose = T){
  
  if(verbose) print("Computing low-dimensional embedding")
  svd_res <- multiomicCCA:::.svd_truncated(mat, K = K, symmetric = F, rescale = F, 
                                           mean_vec = T, sd_vec = F, K_full_rank = F)
  dimred <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
  avg_mat <- .compute_average_mat(dimred, metacell_clustering)
  l2_vec <- apply(avg_mat, 1, multiomicCCA:::.l2norm)
  avg_mat <- multiomicCCA:::.mult_vec_mat(1/l2_vec, avg_mat)
  
  if(verbose) print("Computing distance matrix")
  dist_mat <- as.matrix(stats::dist(avg_mat, method = "euclidean"))
  
  if(verbose) print("Computing bandwidth")
  bandwidth_vec <- .compute_bandwidth(dist_mat, clustering_hierarchy)
  kernel_mat <- matrix(NA, nrow = nrow(dist_mat), ncol = ncol(dist_mat))
  
  if(verbose) print("Computing kernel matrix")
  for(i in 1:nrow(kernel_mat)){
    kernel_mat[i,] <- exp(-dist_mat[i,]/bandwidth_vec[i])
  }
  
  if(verbose) print("Symmeterizing kernel matrix")
  if (symmetrize_func[1] == "max"){
    kernel_mat <- pmax(kernel_mat, t(kernel_mat))
  } else if(symmetrize_func[1] == "min"){
    kernel_mat <- pmin(kernel_mat, t(kernel_mat))
  } else {
    kernel_mat <- (kernel_mat + t(kernel_mat))/2
  }
  
  list(dist_mat = dist_mat, kernel_mat = kernel_mat)
}

compute_min_embedding <- function(kernel_mat_1, 
                                  kernel_mat_2,
                                  K = min(100, nrow(kernel_mat_1)),
                                  entrywise_func = c("max", "min"),
                                  verbose = T){
  stopifnot(all(dim(kernel_mat_1) == dim(kernel_mat_2)))
  n <- nrow(kernel_mat_1)
  
  # entrywise-max
  if(verbose) print("Computing entry-wise maximum")
  if(entrywise_func[1] == "max"){
    kernel_mat <- pmax(kernel_mat_1, kernel_mat_2)
  } else {
    kernel_mat <- pmin(kernel_mat_1, kernel_mat_2)
  }
  
  # projecting to become valid kernel matrix
  if(verbose) print("Projecting into space of PSDs")
  svd_res <- multiomicCCA:::.svd_truncated(kernel_mat, 
                                           K = K, 
                                           symmetric = F,
                                           rescale = F, 
                                           mean_vec = T, 
                                           sd_vec = F, 
                                           K_full_rank = F)
  if(any(svd_res$d <= 0)){
    idx <- which(svd_res$d <= 0)
    svd_res$u <- svd_res$u[,-idx,drop = F]
    svd_res$d <- svd_res$d[-idx]
  }
  kernel_mat <- tcrossprod(multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d), svd_res$u)
  
  # geometric centering
  if(verbose) print("Centering kernel matrix")
  centering_mat <- matrix(-1/n, nrow = n, ncol = n)
  diag(centering_mat) <- 1-1/n
  kernel_mat <- centering_mat %*% kernel_mat %*% centering_mat
  
  # compute embedding
  if(verbose) print("Computing final embedding")
  svd_res <- multiomicCCA:::.svd_truncated(kernel_mat, 
                                           K = K, 
                                           symmetric = F,
                                           rescale = F, 
                                           mean_vec = T, 
                                           sd_vec = F, 
                                           K_full_rank = F)
  dimred <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
}

###############################

.compute_average_mat <- function(mat, metacell_clustering){
  stopifnot(is.list(metacell_clustering))
  
  t(sapply(metacell_clustering, function(vec){
    matrixStats::colSums2(mat[vec,,drop = F])
  }))
}

.compute_bandwidth <- function(dist_mat, clustering_hierarchy){
  tmp <- unlist(clustering_hierarchy); tmp <- tmp[!is.na(tmp)]
  stopifnot(length(unique(tmp)) == length(tmp), 
            all(tmp %% 1 == 0),
            all(sort(tmp) == 1:nrow(dist_mat)),
            nrow(dist_mat) == ncol(dist_mat))
  bandwidth_vec <- rep(NA, max(tmp))
  
  for(i in 1:nrow(dist_mat)){
    k <- which(sapply(clustering_hierarchy, function(vec){i %in% vec}))
    stopifnot(length(k) == 1)
    
    bandwidth_vec[i] <- max(dist_mat[i, clustering_hierarchy[[k]]])
    if(bandwidth_vec[i] == 0){
      bandwidth_vec[i] <- min(dist_mat[i, -clustering_hierarchy[[k]]])/2
    }
  }
  
  stopifnot(all(!is.na(bandwidth_vec)))
  bandwidth_vec
}