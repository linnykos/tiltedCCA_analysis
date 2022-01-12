form_kernel_matrix <- function(mat, K, 
                               metacell_clustering, 
                               clustering_hierarchy,
                               similarity_type = c("euclidean", "kernel")){
  stopifnot(length(similarity_type) == 1, similarity_type == "kernel") # other one is not coded yet
  
  svd_res <- multiomicCCA:::.svd_truncated(mat, K = K, symmetric = F, rescale = F, 
                                           mean_vec = T, sd_vec = F, K_full_rank = F)
  dimred <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
  avg_mat <- .compute_average_mat(dimred, metacell_clustering)
  dist_mat <- as.matrix(stats::dist(avg_mat, method = "euclidean"))
  
  bandwidth_vec <- .compute_bandwidth(dist_mat, clustering_hierarchy)
  kernel_mat <- matrix(NA, nrow = nrow(dist_mat), ncol = ncol(dist_mat))
  
  for(i in 1:nrow(kernel_mat)){
    kernel_mat[i,] <- exp(-dist_mat[i,]/bandwidth_vec[i])
  }
  
  kernel_mat <- (kernel_mat + t(kernel_mat))/2
  kernel_mat
}

compute_min_embedding <- function(kernel_mat_1, 
                                  kernel_mat_2,
                                  K = min(100, nrow(kernel_mat_1))){
  stopifnot(all(dim(kernel_mat_1) == dim(kernel_mat_2)))
  
  kernel_mat <- pmax(kernel_mat_1, kernel_mat_2)
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
  
  # [[note to self: should be consider the laplacian?]]
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
            all(tmp %% 1) == 0,
            length(tmp) == nrow(dist_mat),
            nrow(dist_mat) == ncol(dist_mat))
  bandwidth_vec <- rep(NA, max(tmp))
  
  for(idx in clustering_hierarchy){
    if(all(is.na(idx))) next()
    
    max_val <- max(dist_mat[idx,idx])
    bandwidth_vec[idx] <- max_val
  }
  
  stopifnot(all(!is.na(bandwidth_vec)))
  bandwidth_vec
}
