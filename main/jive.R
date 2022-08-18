## see https://github.com/cran/r.jive/blob/master/R/jive.r
## also based on Section 3 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3671601/bin/NIHMS444793-supplement-Supplemental_Article.pdf

jive <- function(mat_1, mat_2, 
                 common_r, 
                 r_1,
                 r_2,
                 max_iter = 100){
  stopifnot(nrow(mat_1) == nrow(mat_2))
  
  mat_1 <- scale(mat_1, center = T, scale = T)
  mat_2 <- scale(mat_2, center = T, scale = T)
  n <- nrow(mat_1)
  p1 <- ncol(mat_1)
  p2 <- ncol(mat_2)
  
  mat <- cbind(mat_1, mat_2)
  obj_vec <- numeric(0)
  iter <- 1
  
  while(iter <= max_iter){
    print(iter)
    if(iter > 4 && abs(obj_vec[iter-2] - obj_vec[iter-1])/obj_vec[iter-2] <= 1e-6){
      break()
    }
    
    svd_res <- irlba::irlba(mat, nv = common_r)
    embedding <- tiltedCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
    pred_mat <- tcrossprod(embedding, svd_res$v)
    pred_mat_1 <- pred_mat[,1:p1]
    pred_mat_2 <- pred_mat[,(p1+1):(p1+p2)]
    
    resid_1 <- mat_1 - pred_mat_1
    resid_2 <- mat_2 - pred_mat_2
    
    proj_1 <- (diag(n) - tcrossprod(svd_res$u)) %*% resid_1
    proj_2 <- (diag(n) - tcrossprod(svd_res$u)) %*% resid_2
    
    a_svd_1 <- svd_res <- irlba::irlba(proj_1, nv = r_1)
    a_embedding_1 <- tiltedCCA:::.mult_mat_vec(a_svd_1$u, a_svd_1$d)
    a_mat_1 <- tcrossprod(a_embedding_1, a_svd_1$v)
    a_svd_2 <- svd_res <- irlba::irlba(proj_2, nv = r_2)
    a_embedding_2 <- tiltedCCA:::.mult_mat_vec(a_svd_2$u, a_svd_2$d)
    a_mat_2 <- tcrossprod(a_embedding_2, a_svd_2$v)
    
    mat <- cbind(mat_1 - a_mat_1, mat_2 - a_mat_2)
    obj_vec <- c(obj_vec, norm(resid_1, "F")^2 + norm(resid_2, "F")^2)
    iter <- iter + 1
  }
  
  list(embedding = embedding, 
       a_embedding_1 = a_embedding_1,
       a_embedding_2 = a_embedding_2,
       iter = iter,
       obj_vec = obj_vec)
}