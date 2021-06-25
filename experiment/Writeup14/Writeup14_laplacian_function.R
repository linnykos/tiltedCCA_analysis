lap_eig <- function(mat, k_max = 200){
  res <- RSpectra::eigs(mat, k = k_max)
  res$values <- Re(res$values)
  res$vectors <- Re(res$vectors)
  inv_vectors <- MASS::ginv(res$vectors)
  
  approx_mat <- multiomicCCA:::.mult_mat_vec(res$vectors, res$values) %*% inv_vectors
  print(sqrt(sum((mat - approx_mat)^2))/sqrt(sum(mat^2)))
  
  tmp <- multiomicCCA:::.mult_mat_vec(res$vectors[,-1], res$values[-1])
  for(j in 1:ncol(tmp)){
    max_pos <- max(pmax(tmp[,j], 0))
    max_neg <- abs(min(pmin(tmp[,j],0)))
    if(max_neg > max_pos) tmp[,j] <- -tmp[,j]
  }
  
  tmp
}

compute_lap <- function(g, k_max = 100, normalize = T,
                        rowname_vec, colname_vec){
  n <- nrow(g)
  if(normalize){
    deg_vec <- sparseMatrixStats::rowSums2(g)
    invdeg_mat <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1/deg_vec)
    g <- invdeg_mat %*% g %*% invdeg_mat
  }
  
  deg_vec2 <- sparseMatrixStats::rowSums2(g)
  invdeg_mat2 <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1/deg_vec2)
  lap_mat <- invdeg_mat2 %*% g 
  
  basis <- lap_eig(lap_mat, k_max = k_max)
  rownames(basis) <- rowname_vec
  colnames(basis) <- colname_vec
  
  basis
}

compute_smooth_signal <- function(vec, basis){
  lm_res <- stats::lm(vec ~ basis) 
  list(pred_vec = lm_res$fitted.values, r_squared = summary(lm_res)$adj.r.squared)
}
