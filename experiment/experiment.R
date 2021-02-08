rm(list=ls())
library(Seurat)
source("../multiomicCCA_analysis/experiment/Writeup11/Writeup11_simulation_functions.R")

set.seed(10)
B_mat <- matrix(c(0.9, 0.4, 0.1,
                  0.4, 0.9, 0.1,
                  0.1, 0.1, 0.3), 3, 3)
K <- ncol(B_mat); n_clust <- 100; rho <- 1

true_membership_vec <- rep(1:4, each = n_clust)
n <- length(true_membership_vec)

membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, 2*n_clust))
score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)

membership_vec <- c(rep(3, 2*n_clust), rep(1, n_clust), rep(2, n_clust))
score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)

set.seed(10)
p_1 <- 20; p_2 <- 40
coef_mat_1 <- matrix(stats::rnorm(K*p_1), K, p_1)
coef_mat_2 <- matrix(stats::rnorm(K*p_2), K, p_2)

set.seed(10)
dat <- generate_data_dcca(score_1, score_2, coef_mat_1, coef_mat_2)

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

dcca_decomp2 <- dcca_variance_decomposition(dcca_res, rank_c = K, verbose = F)

################################

tol <- 0.01
common_mat <- dcca_decomp2$common_mat_1
mat <- dcca_decomp$common_mat_1+dcca_decomp$distinct_mat_1
tmp <- common_mat/mat
tmp <- pmin(pmax(tmp, tol), 1)
quantile(tmp)
svd_res <- .svd_truncated(tmp, K = 1)
if(sum(svd_res$u<0) + sum(svd_res$v<0) > sum(svd_res$u>0) + sum(svd_res$v>0)){
  svd_res$u <- -svd_res$u; svd_res$v <- -svd_res$v
}
alpha_vec <- as.numeric(svd_res$u*sqrt(svd_res$d))
quantile(alpha_vec); alpha_vec <- pmin(pmax(alpha_vec, tol), 1); quantile(alpha_vec)
beta_vec <- as.numeric(svd_res$v*sqrt(svd_res$d))
quantile(beta_vec); beta_vec <- pmin(pmax(beta_vec, tol), 1); quantile(beta_vec)

# now alternating minimization
n <- length(alpha_vec); p <- length(beta_vec)
obj_vec <- rep(NA, iter_max)
iter_max <- 100
for(iter in 1:iter_max){
  print(iter)
  obj_vec[iter] <- .l2norm(common_mat - .mult_vec_mat(alpha_vec, .mult_mat_vec(mat, beta_vec)))^2
  
  tmp <- .mult_mat_vec(mat, beta_vec)
  alpha_vec <- sapply(1:n, function(i){
    as.numeric(common_mat[i,] %*% tmp[i,])/crossprod(tmp[i,])
  })
  alpha_vec <- pmin(pmax(alpha_vec, tol), 1)
  
  tmp <- .mult_vec_mat(alpha_vec, mat)
  beta_vec <- sapply(1:p, function(j){
    as.numeric(common_mat[,j] %*% tmp[,j])/crossprod(tmp[,j])
  })
  beta_vec <- pmin(pmax(beta_vec, tol), 1)
}
plot(obj_vec)
quantile(alpha_vec)
quantile(beta_vec)
plot(alpha_vec)

obj_vec[iter_max]/.l2norm(mat)^2
