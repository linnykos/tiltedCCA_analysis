form_dist_matrix <- function(mat, K, 
                             metacell_clustering, 
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
  dist_mat^2
}

compute_min_embedding <- function(dist_mat_1, 
                                  dist_mat_2,
                                  sing_val_cutoff = 0.01,
                                  verbose = T){
  stopifnot(all(dim(kernel_mat_1) == dim(kernel_mat_2)))
  n <- nrow(kernel_mat_1)
  
  # entrywise-max
  if(verbose) print("Computing entry-wise minimum")
  dist_mat <- pmin(dist_mat_1, dist_mat_2)
  
  # projecting to become valid kernel matrix
  if(verbose) print("Projecting into space of distance matries")
  dist_mat <- .projection_to_edm_cone(dist_mat)
  
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
.projection_to_edm_cone <- function(mat, iter_max = 100, tol = 1e-4){
  stopifnot(is.matrix(mat), nrow(mat) == ncol(mat))
  
  n <- nrow(mat)
  centering_mat <- matrix(-1/n, nrow = n, ncol = n)
  diag(centering_mat) <- 1-1/n
  
  iter <- 1
  prev_mat <- mat
  
  while(TRUE){
    if(iter > 1 && sum(abs(prev_mat - mat)) <= tol) break()
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
  }
  
  mat
}
