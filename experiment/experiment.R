rm(list=ls())
library(Seurat)
source("../multiomicCCA_analysis/experiment/Writeup12/Writeup12_simulation_functions.R")

set.seed(10)
n_clust <- 100
B_mat <- matrix(c(0.9, 0.2, 0.1, 
                  0.2, 0.9, 0.1,
                  0.1, 0.1, 0.5), 3, 3)
K <- ncol(B_mat)
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
n <- length(membership_vec); true_membership_vec <- membership_vec
score_1 <- generate_sbm_orthogonal(B_mat, membership_vec)
score_2 <- generate_sbm_orthogonal(B_mat, membership_vec)

set.seed(10)
p_1 <- 20; p_2 <- 40
coef_mat_1 <- matrix(stats::rnorm(K*p_1), K, p_1)
coef_mat_2 <- matrix(stats::rnorm(K*p_2), K, p_2)

set.seed(10)
dat <- generate_data_dcca(score_1, score_2, coef_mat_1, coef_mat_2)

# set.seed(10)
# dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
# dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

###############################

mat_1 = dat$mat_1; mat_2 = dat$mat_2;
rank_1 = K; rank_2 = K; meta_clustering = NA
num_neigh = max(round(nrow(mat_1)/20), 40)
apply_shrinkage = F
reorthogonalize = F
verbose = T

stopifnot(nrow(mat_1) == nrow(mat_2), 
          rank_1 <= min(dim(mat_1)), rank_2 <= min(dim(mat_2)))
n <- nrow(mat_1)
if(verbose) print("D-CCA: Rescaling matrices")
mat_1 <- scale(mat_1, center = T, scale = F)
mat_2 <- scale(mat_2, center = T, scale = F)

if(verbose) print(paste0("D-CCA: Starting matrix shrinkage"))
if(apply_shrinkage) svd_1 <- .spoet(mat_1, rank_1) else svd_1 <- .svd_truncated(mat_1, rank_1)
if(apply_shrinkage) svd_2 <- .spoet(mat_2, rank_2) else svd_2 <- .svd_truncated(mat_2, rank_2)

svd_1 <- .check_svd(svd_1); svd_2 <- .check_svd(svd_2)

msg <- " (all cells)"
if(verbose) print(paste0("D-CCA", msg, ": Computing CCA"))
cca_res <- .cca(svd_1, svd_2)




