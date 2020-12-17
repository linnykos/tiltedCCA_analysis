rm(list=ls())
set.seed(5)
n <- 100; K <- 2
common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)

p1 <- 5; p2 <- 10
transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)

mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = 0.01*diag(p1)), center = T, scale = F)
mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = 0.01*diag(p2)), center = T, scale = F)

# dcca_res <- dcca(mat_1, mat_2, rank_1 = K, rank_2 = K, rank_12 = K)

######
rank_1 = K; rank_2 = K; rank_12 = K; enforce_rank = T

stopifnot(nrow(mat_1) == nrow(mat_2), rank_12 <= min(c(rank_1, rank_2)), 
          rank_1 <= min(dim(mat_1)), rank_2 <= min(dim(mat_2)))
n <- nrow(mat_1)
mat_1 <- scale(mat_1, center = T, scale = F)
mat_2 <- scale(mat_2, center = T, scale = F)

if(enforce_rank | nrow(mat_1) < 2*ncol(mat_1)) mat_1 <- .spoet(mat_1, rank_1)
if(enforce_rank | nrow(mat_2) < 2*ncol(mat_2)) mat_2 <- .spoet(mat_2, rank_2)

cov_1 <- stats::cov(mat_1) * (n-1)/n
cov_2 <- stats::cov(mat_2) * (n-1)/n
cov_12 <- crossprod(mat_1, mat_2)/n
full_rank <- Matrix::rankMatrix(cov_12)

cca_res <- .cca(cov_1, cov_2, cov_12, K = K)

score_1 <- mat_1 %*% cca_res$factor_1/sqrt(n)
score_2 <- mat_2 %*% cca_res$factor_2/sqrt(n)

crossprod(score_1)
crossprod(score_2)


cca_res <- .cca(cov_1, cov_2, cov_12, K = full_rank)
score_1 <- mat_1 %*% cca_res$factor_1 #note: this is unnormalized
score_2 <- mat_2 %*% cca_res$factor_2

R_vec <- sapply(cca_res$obj_vec, function(x){1-sqrt((1-x)/(1+x))})
