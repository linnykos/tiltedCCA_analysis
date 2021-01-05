rm(list=ls())
library(Seurat)
source("../multiomicCCA_analysis/experiment/Writeup10/Writeup10_simulation_functions.R")

# try a setting where the common loading separates (12) from (34), 
# distinct 1 separates (1) from (2) [with high weights] and 
# distinct 2 seperates (3) from (4) analogously
# common space will have very small weight
set.seed(10)
B_mat <- matrix(c(0.9, 0.4, 
                  0.4, 0.9), 2, 2)
K <- ncol(B_mat); n_clust <- 100; n <- 5*n_clust; rho <- 0.3

true_membership_vec <- rep(1:4, each = n_clust)
membership_vec <- c(rep(1, 2*n_clust), rep(2, 2*n_clust))
score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)

membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(1, n_clust), rep(2, n_clust))
score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)

set.seed(10)
p_1 <- 20; p_2 <- 40
coef_mat_1 <- matrix(stats::rnorm(K*p_1), K, p_1)
coef_mat_2 <- matrix(stats::rnorm(K*p_2), K, p_2)

set.seed(10)
dat <- generate_data(score_1, score_2, coef_mat_1, coef_mat_2)

par(mfrow = c(1,3))
image(t(dat$common_score))
image(t(dat$distinct_score_1))
image(t(dat$distinct_score_2))

