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
dat <- generate_data_dcca(score_1, score_2, coef_mat_1, coef_mat_2, noise_func = function(x){x}); class(dat) <- "dcca_decomp"

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F, reorthogonalize = T)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

png("../../out/simulation/Writeup11/Writeup11_simulation2_reorth_true_score.png", height = 1200, width = 2500, units = "px", res = 300)
plot_scores(dcca_res, mode = 1)
graphics.off()

png("../../out/simulation/Writeup11/Writeup11_simulation2_reorth_true_embedding.png", height = 800, width = 2000, units = "px", res = 300)
plot_embeddings(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()

