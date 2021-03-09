rm(list=ls())
library(Seurat)
source("../multiomicCCA_analysis/experiment/Writeup12/Writeup12_simulation_functions.R")

set.seed(10)
n_clust <- 100
B_mat1 <- matrix(c(1, 0, 0, 
                   0, 1, 0,
                   0, 0, 1), 3, 3)
B_mat2 <- matrix(c(1, 1, 0, 
                   1, 1, 0,
                   0, 0, 1), 3, 3)
K <- ncol(B_mat1)
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
n <- length(membership_vec); true_membership_vec <- membership_vec
score_1 <- generate_sbm_orthogonal(B_mat1, membership_vec)
score_2 <- generate_sbm_orthogonal(B_mat2, membership_vec)

set.seed(10)
p_1 <- 20; p_2 <- 40
coef_mat_1 <- matrix(stats::rnorm(K*p_1), K, p_1)
coef_mat_2 <- matrix(stats::rnorm(K*p_2), K, p_2)

set.seed(10)
dat <- generate_data_dcca(score_1, score_2, coef_mat_1, coef_mat_2, 
                          noise_func = function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat, sd = 2.5), nrow(mat), ncol(mat))})

mat_1 = dat$mat_1; mat_2 = dat$mat_2
rank_1 = K; rank_2 = K
apply_shrinkage = F
verbose = F
n <- nrow(mat_1)
if(verbose) print("D-CCA: Rescaling matrices")
mat_1 <- scale(mat_1, center = T, scale = F)
mat_2 <- scale(mat_2, center = T, scale = F)

if(verbose) print(paste0("D-CCA: Starting matrix shrinkage"))
if(apply_shrinkage) svd_1 <- .spoet(mat_1, rank_1) else svd_1 <- .svd_truncated(mat_1, rank_1)
if(apply_shrinkage) svd_2 <- .spoet(mat_2, rank_2) else svd_2 <- .svd_truncated(mat_2, rank_2)

svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)
cca_res <- .cca(svd_1, svd_2)

num_neigh = max(round(nrow(svd_1$u)/20), 40)
check_alignment = T; reorthogonalize = F; verbose = T
full_rank <- length(cca_res$obj_vec)
n <- nrow(svd_1$u)

tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
score_1 <- tmp$score_1; score_2 <- tmp$score_2
stopifnot(ncol(score_1) == length(svd_1$d), ncol(score_2) == length(svd_2$d),
          nrow(score_1) == nrow(score_2))

tmp <- .cca(score_1, score_2, rank_1 = ncol(score_1), rank_2 = ncol(score_2), return_scores = T)
score_1 <- tmp$score_1; score_2 <- tmp$score_2
stopifnot(is.matrix(score_1), is.matrix(score_2))
obj_vec <- diag(crossprod(score_1, score_2))/n

nn_1 <- RANN::nn2(tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_1$v), k = 50)
nn_2 <- RANN::nn2(tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_2$v), k = 50)

################

rank_c <- min(ncol(score_1), ncol(score_2))
k = 1
common_perc <- .latent_common_perc(score_1[,k], score_2[,k], nn_1, nn_2)

basis_res <- .representation_2d(score_1[,k], score_2[,k])
vec1 <- basis_res$rep1; vec2 <- basis_res$rep2
ang <- .angle_between_vectors(vec1, vec2)
vec_right <- .rightmost_vector(vec1, vec2)
if(sum(abs(vec_right$vec_right - vec1)) <= sum(abs(vec_right$vec_right - vec1))){
  common_vec <- .angle_from_vector(vec_right$vec_right, (1-common_perc)*ang)
} else {
  common_vec <- .angle_from_vector(vec_right$vec_right, common_perc*ang)
}

a <- .l2norm(common_vec)^2; b <- -as.numeric(common_vec %*% (vec1 + vec2)); c <- as.numeric(vec1 %*% vec2)
common_vec <- .quadratic(a,b,c) * common_vec

(vec1 - common_vec) %*% (vec2 - common_vec)
z1 <- .angle_between_vectors(vec1, common_vec)
z2 <- .angle_between_vectors(vec2, common_vec)
z1; z2

.decomposition_2d(vec1, vec2, common_vec)

