form_dist_matrix <- function(mat, K, 
                             metacell_clustering, 
                             clustering_hierarchy = NA,
                             regularization_quantile = 0.1,
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
  dist_mat <- dist_mat^2
  rownames(dist_mat) <- names(metacell_clustering)
  colnames(dist_mat) <- names(metacell_clustering)
  
  if(!all(is.na(clustering_hierarchy))){
    if(verbose) print("Apply large-clustering regularization")
    
    metacell_affiliation <- sapply(names(metacell_clustering), function(string){
      paste0("large_", strsplit(string, split = "_")[[1]][2])
    })
    
    for(i in 1:length(clustering_hierarchy)){
      if(length(clustering_hierarchy[[i]]) == 1) next()
      combn_mat <- utils::combn(length(clustering_hierarchy[[i]]), 2)
      for(j in 1:ncol(combn_mat)){
        large_1 <- clustering_hierarchy[[i]][combn_mat[1,j]]
        large_2 <- clustering_hierarchy[[i]][combn_mat[2,j]]
        
        idx_1 <- which(metacell_affiliation == large_1)
        idx_2 <- which(metacell_affiliation == large_2)
        
        max_val <- stats::quantile(dist_mat[idx_1, idx_2], probs = regularization_quantile)
        dist_mat[idx_1, idx_2] <- pmin(dist_mat[idx_1, idx_2], max_val)
        dist_mat[idx_2, idx_1] <- pmin(dist_mat[idx_2, idx_1], max_val)
      }
    }
  }
  
  dist_mat
}

compute_min_embedding <- function(dist_mat_1, 
                                  dist_mat_2,
                                  sing_val_cutoff = 0.01,
                                  verbose = T){
  stopifnot(all(dim(dist_mat_1) == dim(dist_mat_2)))
  n <- nrow(dist_mat_1)
  
  # entrywise-max
  if(verbose) print("Computing entry-wise minimum")
  dist_mat <- pmin(dist_mat_1, dist_mat_2)
  
  # projecting to become valid kernel matrix
  if(verbose) print("Projecting into space of distance matries")
  dist_mat <- .projection_to_edm_cone(dist_mat, 
                                      verbose = verbose)
  
  # geometric centering
  if(verbose) print("Convering into gram matrix")
  centering_mat <- matrix(-1/n, nrow = n, ncol = n)
  diag(centering_mat) <- 1-1/n
  gram_mat <- -.5*centering_mat %*% dist_mat %*% centering_mat
  
  # compute embedding
  if(verbose) print("Computing final embedding")
  svd_res <- multiomicCCA:::.svd_truncated(gram_mat, 
                                           K = min(100, nrow(gram_mat)), 
                                           symmetric = F,
                                           rescale = F, 
                                           mean_vec = F, 
                                           sd_vec = F, 
                                           K_full_rank = F)
  
  spectrum_ratio <- abs(diff(svd_res$d))/svd_res$d[-1]
  K <- which(spectrum_ratio <= sing_val_cutoff)[1]
  if(length(K) == 0) {
    K <- length(svd_res$d)
  } 
  dimred <- multiomicCCA:::.mult_mat_vec(svd_res$u[,1:K,drop=F], svd_res$d[1:K])
  
  list(dist_mat = dist_mat, dimred = dimred)
}

###############################

.compute_average_mat <- function(mat, metacell_clustering){
  stopifnot(is.list(metacell_clustering))
  
  t(sapply(metacell_clustering, function(vec){
    matrixStats::colSums2(mat[vec,,drop = F])
  }))
}

# see "Equality relating Euclidean distance cone to positivie definite cone", Dattoro 2008
.projection_to_edm_cone <- function(mat, 
                                    iter_max = 200, 
                                    tol = 1e-4,
                                    verbose = F){
  stopifnot(is.matrix(mat), nrow(mat) == ncol(mat))
  
  n <- nrow(mat)
  centering_mat <- matrix(-1/n, nrow = n, ncol = n)
  diag(centering_mat) <- 1-1/n
  
  iter <- 1
  prev_mat <- mat
  
  while(TRUE){
    if(iter > 1 && sum(abs(prev_mat - mat)) <= tol) break()
    if(iter > iter_max) break()
    if(verbose) print(paste0("Iteration ", iter, " with difference ", round(sum(abs(prev_mat - mat)))))
    prev_mat <- mat
    
    # project to symmetric hollow subspace
    mat <- (mat + t(mat))/2
    diag(mat) <- 0
    
    # project so that negative centered matrix is PSD
    centered_mat <- centering_mat %*% mat %*% centering_mat
    eigen_res <- eigen(centered_mat)
    eigen_res$values[eigen_res$values < 0] <- 0
    centered_mat <- tcrossprod(multiomicCCA:::.mult_mat_vec(eigen_res$vectors, eigen_res$values), eigen_res$vectors)
    
    mat <- mat - centered_mat
    iter <- iter + 1
  }
  
  mat
}
