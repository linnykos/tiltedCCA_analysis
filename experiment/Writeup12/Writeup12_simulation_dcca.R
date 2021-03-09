# easy setting: low distinct 1, low distinct 2 with 3 clusters
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
dat <- generate_data(score_1, score_2, coef_mat_1, coef_mat_2)

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

png("../../out/simulation/Writeup12/Writeup12_simulation1_embedding.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec, observed = F)
graphics.off()

par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F)
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T)

#######################################

# high distinct 1, low distinct 2 with 3 clusters
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
dat <- generate_data(score_1, score_2, coef_mat_1, coef_mat_2, 
                          noise_func = function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat, sd = 2.5), nrow(mat), ncol(mat))})
set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

# plot_embeddings(dcca_decomp, membership_vec = true_membership_vec)
png("../../out/simulation/Writeup12/Writeup12_simulation2_embedding.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec, observed = F)
graphics.off()

par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F)
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T)

##################################

# high distinct 1, high distinct 2 with 3 clusters
set.seed(10)
n_clust <- 100
high <- 0.9; low <- 0.05
B_mat1 <- matrix(c(high, high, low, 
                   high, high, low,
                   low, low, high), 3, 3)
B_mat2 <- matrix(c(high, low, low, 
                   low, high, high,
                   low, high, high), 3, 3)
K <- ncol(B_mat1)
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
n <- length(membership_vec); true_membership_vec <- membership_vec
score_1 <- generate_sbm_orthogonal(B_mat1, membership_vec)
score_2 <- generate_sbm_orthogonal(B_mat2, membership_vec)
#score_2 <- score_1[c((2*n_clust+1):(3*n_clust),1:(2*n_clust)),]

set.seed(10)
p_1 <- 20; p_2 <- 20
coef_mat1 <- matrix(stats::rnorm(K*p_1), K, p_1)
coef_mat2 <- matrix(stats::rnorm(K*p_1), K, p_1)

set.seed(10)
dat <- generate_data(score_1, score_2, coef_mat1, coef_mat2, 
                          noise_func = function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat, sd = 1.5), nrow(mat), ncol(mat))})
set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

png("../../out/simulation/Writeup12/Writeup12_simulation3_embedding.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()

plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec)
plot_scores(dcca_decomp, membership_vec = true_membership_vec)

par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F)
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T)
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T)

####################################3

# high distinct 1, high distinct 2 with 4 clusters. 
# distinct 1 separates (1) from (2) and 
# distinct 2 seperates (3) from (4) analogously

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
dat <- generate_data(score_1, score_2, coef_mat1, coef_mat2, 
                          noise_func = function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat, sd = 0.5), nrow(mat), ncol(mat))})

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec)

par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F)
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T)
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T)

##################################

# high distinct 1, high distinct 2 with 5 clusters
# the common loading separates (1234) from (5)

set.seed(10)
B_mat <- matrix(c(0.9, 0.4, 0.1, 
                  0.4, 0.9, 0.1, 
                  0.1, 0.1, 0.9), 3, 3)
K <- ncol(B_mat); n_clust <- 100; rho <- 1
true_membership_vec <- rep(1:5, each = n_clust)
n <- length(true_membership_vec)

membership_vec <-  c(rep(1, 2*n_clust), rep(2, 2*n_clust), rep(3, n_clust))
score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)

membership_vec <-  c(rep(1, n_clust), rep(2, n_clust), rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)

set.seed(10)
p_1 <- 20; p_2 <- 40
coef_mat_1 <- matrix(stats::rnorm(K*p_1), K, p_1)
coef_mat_2 <- matrix(stats::rnorm(K*p_2), K, p_2)

set.seed(10)
dat <- generate_data(score_1, score_2, coef_mat1, coef_mat2, 
                          noise_func = function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat, sd = 0.5), nrow(mat), ncol(mat))})

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec)

par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F)
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T)
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T)



