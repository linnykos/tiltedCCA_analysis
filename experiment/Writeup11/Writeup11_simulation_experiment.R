rm(list=ls())
library(Seurat)
source("../multiomicCCA_analysis/experiment/Writeup11/Writeup11_simulation_functions.R")

# try a setting where the common loading separates (12) from (34), 
# distinct 1 separates (1) from (2) [with high weights] and 
# distinct 2 seperates (3) from (4) analogously
# common space will have very small weight
set.seed(10)
B_mat <- matrix(c(0.9, 0.4, 0.1,
                  0.4, 0.9, 0.1,
                  0.1, 0.1, 0.5), 3, 3)
K <- ncol(B_mat); n_clust <- 100; 

true_membership_vec <- rep(1:4, each = n_clust)
n <- length(true_membership_vec)
rho <- 1

membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, 2*n_clust))
prob_mat1 <- .compute_prob_mat(rho*B_mat, membership_vec)
membership_vec <- c(rep(3, 2*n_clust), rep(1, n_clust), rep(2, n_clust))
prob_mat2 <- .compute_prob_mat(rho*B_mat, membership_vec)

set.seed(10)
p_1 <- 20; p_2 <- 40
coef_mat_1 <- matrix(stats::rnorm(K*p_1), K, p_1)
coef_mat_2 <- matrix(stats::rnorm(K*p_2), K, p_2)
coef_mat_1 <- apply(coef_mat_1, 2, sort)
coef_mat_2 <- apply(coef_mat_2, 2, sort)
coef_mat_1 <- coef_mat_1/max(svd(coef_mat_1)$d[1])*0.9
coef_mat_2 <- coef_mat_2/max(svd(coef_mat_2)$d[1])*0.9

########

svd_1 <- .svd_truncated(prob_mat1, K)
svd_2 <- .svd_truncated(prob_mat2, K)
# zz <- t(svd_1$u) %*% svd_2$u; Matrix::rankMatrix(zz)
# image(zz)

################

mat_1 <- svd_1$u %*% diag(svd_1$d) %*% coef_mat_1
mat_2 <- svd_2$u %*% diag(svd_2$d) %*% coef_mat_2

# par(mfrow=c(1,2))
# image(t(mat_1)); image(t(mat_2))

svd_full_1 <- .svd_truncated(mat_1, K)
svd_full_2 <- .svd_truncated(mat_2, K)

# par(mfrow=c(1,2))
# image(t(svd_full_1$u)); image(t(svd_full_2$u))

####################

set.seed(10)
dmca_res <- dmca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
res <- dmca_decomposition(dmca_res, verbose = F)
plot_scores(res, mode = 2)
plot_embeddings(res, true_membership_vec)

# why is there signal in the common space...?
par(mfrow=c(1,2))
image(t(res$common_score_1 + res$distinct_score_1))
image(t(res$common_score_2 + res$distinct_score_2))
image(crossprod(mat_1, mat_2))



set.seed(10)
dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
res <- dcca_decomposition(dcca_res, rank_c = 2, verbose = F)
plot_scores(res, mode = 1)
plot_embeddings(res, true_membership_vec)

par(mfrow=c(1,2))
image(t(res$common_score + res$distinct_score_1))
image(t(res$common_score + res$distinct_score_2))




