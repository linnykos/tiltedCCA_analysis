rm(list=ls())
set.seed(5)
n <- 100; K <- 3
common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)

p1 <- 5; p2 <- 10
transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)

mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)

# dcca_res <- dcca(mat_1, mat_2, rank_1 = K, rank_2 = K, rank_12 = K)

######
rank_1 = K; rank_2 = K; rank_12 = K; enforce_rank = T; verbose = T

n <- nrow(mat_1)
mat_1 <- scale(mat_1, center = T, scale = F)
mat_2 <- scale(mat_2, center = T, scale = F)

if(verbose) print("D-CCA: Starting matrix shrinkage")
if(enforce_rank | nrow(mat_1) < 2*ncol(mat_1)) mat_1 <- .spoet(mat_1, rank_1)
if(enforce_rank | nrow(mat_2) < 2*ncol(mat_2)) mat_2 <- .spoet(mat_2, rank_2)

###########

svd_res_1 <- .svd_truncated(mat_1, rank_1)
svd_res_2 <- .svd_truncated(mat_2, rank_2)

cov_1 <- stats::cov(mat_1) * (n-1)/n
cov_2 <- stats::cov(mat_2) * (n-1)/n
cov_12 <- crossprod(mat_1, mat_2)/n
full_rank <- Matrix::rankMatrix(cov_12)

cov_1_invhalf <- .mult_mat_vec(svd_res_1$v, 1/svd_res_1$d)
cov_2_invhalf <- .mult_mat_vec(svd_res_2$v, 1/svd_res_2$d)

# test that it really is an inverse
sum(abs(diag(K) - n*crossprod(cov_1_invhalf, cov_1 %*% cov_1_invhalf)))
sum(abs(diag(K) - n*crossprod(cov_2_invhalf, cov_2 %*% cov_2_invhalf)))

prod_v1 <- n*crossprod(cov_1_invhalf, cov_12) %*% cov_2_invhalf
prod_v2 <- crossprod(svd_res_1$u, svd_res_2$u)
sum(abs(prod_v1 - prod_v2))

