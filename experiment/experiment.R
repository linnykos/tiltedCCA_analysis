rm(list=ls())
x = 1
n <- 200; K <- 2
common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)

p1 <- 5; p2 <- 10
transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)

mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)
nc <- 20
meta_clustering <- stats::kmeans(mat_1, centers = nc)$cluster

mat_1 <- scale(mat_1, center = T, scale = F)
mat_2 <- scale(mat_2, center = T, scale = F)

svd_1 <- .svd_truncated(mat_1, K); svd_2 <- .svd_truncated(mat_2, K)
svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)

num_meta <- max(meta_clustering)

mat_1_meta <- t(sapply(1:num_meta, function(x){
  idx <- which(meta_clustering == x)
  apply(mat_1[idx,,drop = F], 2, mean)
}))

mat_2_meta <- t(sapply(1:num_meta, function(x){
  idx <- which(meta_clustering == x)
  apply(mat_2[idx,,drop = F], 2, mean)
}))

cca_res <- .cca(mat_1_meta, mat_2_meta, rank_1 = K, rank_2 = K)

res <- .dcca_common_factors(svd_1, svd_2, cca_res, check_alignment = T, verbose = F)

## 

full_rank <- length(cca_res$obj_vec)
tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
score_1 <- tmp$score_1; score_2 <- tmp$score_2
tmp <- .reparameterize(score_1, score_2)
score_1 <- tmp$mat_1; score_2 <- tmp$mat_2

crossprod(score_1)
crossprod(score_2)
crossprod(score_1, score_2)
diag(crossprod(score_1, score_2))/n
res$cca_obj

residual_1 <- score_1 - res$common_factors
residual_2 <- score_2 - res$common_factors

prod_mat <- t(residual_1) %*% residual_2
sum(abs(prod_mat)) <= 1e-6
