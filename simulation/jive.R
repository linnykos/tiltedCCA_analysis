jive <- function(mat_1, mat_2, r, max_iter = 10){
  mat_1 <- scale(mat_1, center = T, scale = T)
  mat_2 <- scale(mat_2, center = T, scale = T)
  n <- nrow(mat_1)
  p1 <- ncol(mat_1)
  p2 <- ncol(mat_2)
  
  mat <- cbind(mat_1, mat_2)
  embedding_old <- numeric(0)
  iter <- 1
  
  while(TRUE){
    if(is.matrix(embedding_old) && iter > 4){
      if(svd(embedding_old - embedding)$d[1] <= 1e-3) break()
    } 
    if(iter >= max_iter) break()
    
    svd_res <- irlba::irlba(mat, nv = r)
    embedding <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
    pred_mat <- tcrossprod(embedding, svd_res$v)
    pred_mat_1 <- pred_mat[,1:p1]
    pred_mat_2 <- pred_mat[,(p1+1):(p1+p2)]
    
    resid_1 <- mat_1 - pred_mat_1
    resid_2 <- mat_2 - pred_mat_2
    
    proj_1 <- (diag(n) - tcrossprod(svd_res$u)) %*% resid_1
    proj_2 <- (diag(n) - tcrossprod(svd_res$u)) %*% resid_2
    
    mat <- cbind(mat_1 - proj_1, mat_2 - proj_2)
    embedding_old <- embedding
    iter <- iter + 1
  }

  list(embedding = embedding, iter = iter)
}