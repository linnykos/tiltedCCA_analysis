zongming_embedding <- function(mat_1, mat_2,
                               dims_1, dims_2, 
                               nn){
  rank_1 <- max(dims_1); rank_2 <- max(dims_2)
  n <- nrow(mat_1)
  
  svd_1 <- multiomicCCA:::.svd_truncated(mat_1, K = rank_1, symmetric = F, rescale = F, 
                          mean_vec = F, sd_vec = F, K_full_rank = F)
  svd_2 <- multiomicCCA:::.svd_truncated(mat_2, K = rank_2, symmetric = F, rescale = F, 
                          mean_vec = F, sd_vec = F, K_full_rank = F)
  
  svd_1 <- multiomicCCA:::.check_svd(svd_1, dims = dims_1)
  svd_2 <- multiomicCCA:::.check_svd(svd_2, dims = dims_2)
  
  dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)
  dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d)
  
  print("Computing NNs")
  nn_res_1 <- RANN::nn2(dimred_1, k = nn+1)$nn.idx[,-1]
  nn_res_2 <- RANN::nn2(dimred_2, k = nn+1)$nn.idx[,-1]
  
  print("Joining NNs")
  nn_list <- lapply(1:n, function(i){
    if(i %% floor(n/10) == 0) cat('*')
    sort(unique(c(nn_res_1[i,], nn_res_2[i,])))
  })
  
  print("Computing kernels for Modality 1")
  kernel_list_1 <- lapply(1:n, function(i){
    if(i %% floor(n/10) == 0) cat('*')
    dist_vec <- sapply(nn_list[[i]], function(j){
      multiomicCCA:::.l2norm(dimred_1[i,] - dimred_1[j,])
    })
    min_val <- min(dist_vec)
    max_val <- max(dist_vec)
    
    exp(-(dist_vec - min_val)/max_val)
  })
  
  print("Computing kernels for Modality 2")
  kernel_list_2 <- lapply(1:n, function(i){
    if(i %% floor(n/10) == 0) cat('*')
    dist_vec <- sapply(nn_list[[i]], function(j){
      multiomicCCA:::.l2norm(dimred_2[i,] - dimred_2[j,])
    })
    min_val <- min(dist_vec)
    max_val <- max(dist_vec)
    
    exp(-(dist_vec - min_val)/max_val)
  })
  
  print("Computing min kernel")
  dist_list <- lapply(1:n, function(i){
    if(i %% floor(n/10) == 0) cat('*')
    pmin(kernel_list_1[[i]], kernel_list_2[[i]])
  })
  
  print("Forming matrix")
  i_vec <- rep(1:n, times = sapply(nn_list, length))
  j_vec <- unlist(nn_list)
  x_vec <- unlist(dist_list)
  
  sparse_mat <- Matrix::sparseMatrix(i = i_vec,
                                     j = j_vec,
                                     x = x_vec,
                                     dims = c(n,n),
                                     repr = "C")
  colnames(sparse_mat) <- rownames(mat_1)
  rownames(sparse_mat) <- rownames(mat_1)
  
  sparse_mat
}