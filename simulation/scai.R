scai <- function(mat_1, mat_2, r, 
                 gamma = 0, max_iter = 100){
 
  n <- nrow(mat_1)
  p1 <- ncol(mat_1)
  p2 <- ncol(mat_2)
  
  # initialization
  X1t <- t(mat_1); X2t <- t(mat_2)
  H <- abs(matrix(stats::runif(n * r), n, r))
  W1 <- abs(matrix(stats::runif(p1 * r), p1, r))
  W2 <- abs(matrix(stats::runif(p2 * r), p2, r))
  obj_vec <- numeric(0)
  iter <- 1
  
  while(iter <= max_iter){
    print(iter)
    if(iter > 4 && abs(obj_vec[iter-2] - obj_vec[iter-1])/obj_vec[iter-2] <= 1e-6){
      break()
    }
    
    W1 <- W1 * (X1t %*% H)/(W1 %*% crossprod(H))
    W2 <- W2 * (X2t %*% H)/(W2 %*% crossprod(H))
    
    num <- mat_1 %*% W1 + mat_2 %*% W2
    denom <- H %*% (crossprod(W1) + crossprod(W2) + gamma * matrix(1, r, r))
    H <- H * num/denom
    
    obj_val <- norm(mat_1 - tcrossprod(H, W1), "F")^2
    obj_val <- obj_val + norm(mat_2 - tcrossprod(H, W2), "F")^2
    obj_val <- obj_val + gamma * sum(matrixStats::rowSums2(H)^2)
    obj_vec <- c(obj_vec, obj_val)
    iter <- iter + 1
  }
  
  list(W1 = W1, W2 = W2, H = H,
       obj_vec = obj_vec)
}