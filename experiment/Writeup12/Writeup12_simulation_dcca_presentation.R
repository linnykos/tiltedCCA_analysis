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

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

png("../../out/simulation/Writeup12/Writeup12_simulation1_embedding.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()

par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, membership_vec = true_membership_vec)

#######################################

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
set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
# dcca_decomp2 <- dcca_variance_decomposition(dcca_res, rank_c = K, verbose = F)

# plot_embeddings(dcca_decomp, membership_vec = true_membership_vec)
png("../../out/simulation/Writeup12/Writeup12_simulation2_embedding.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()

par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, membership_vec = true_membership_vec)


# plot_scores(dcca_decomp)

# res <- custom_explained_variance(dcca_decomp)
# weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
# par(mfrow = c(1,1), mar = c(4, 4, 3.5, 0.5))
# information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F,
#                  main = "Simulation 2 (L2 Ratio)")
# 
# res <- custom_explained_variance2(dcca_decomp, true_membership_vec)
# weight_mat <- data.frame(Mode_1 = res$weight_vec_1, Mode_2 = res$weight_vec_2, cell_type = 1:length(unique(true_membership_vec)))
# par(mfrow = c(1,1), mar = c(4, 4, 3.5, 0.5))
# information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F,
#                  main = "Simulation 2 (R-Squared)")
# 
# res <- explained_variance(dcca_decomp2, verbose = F)
# weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
# par(mfrow = c(1,1), mar = c(4, 4, 3.5, 0.5))
# information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F,
#                  main = "Simulation 2 (Projection)")
# title(xlab="Cell-type", mgp=c(1,1,0))

##################################


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
dat <- generate_data_dcca(score_1, score_2, coef_mat1, coef_mat2, 
                          noise_func = function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat, sd = 1.5), nrow(mat), ncol(mat))})
set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

png("../../out/simulation/Writeup12/Writeup12_simulation3_embedding.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()

par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, membership_vec = true_membership_vec)

par(mfrow = c(1,3))
image(t(dcca_res$common_score))
image(t(dcca_res$distinct_score_1))
image(t(dcca_res$distinct_score_2))

set.seed(10)
svd_list <- extract_svd_embedding(dcca_decomp); noise_val = 0.05
par(mfrow = c(1,3))
set.seed(10)
tmp <- extract_umap_embedding(svd_list, common_1 = T, common_2 = F, distinct_1 = F, distinct_2 = F, 
                              vis_param = dcca_decomp$vis_param, noise_val = noise_val)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Common view",
     xlab = "UMAP 1", ylab = "UMAP 2")
set.seed(10)
tmp <- extract_umap_embedding(svd_list, common_1 = F, common_2 = F, distinct_1 = T, distinct_2 = F, 
                              vis_param = dcca_decomp$vis_param, noise_val = noise_val)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Distinct view",
     xlab = "UMAP 1", ylab = "UMAP 2")
set.seed(10)
tmp <- extract_umap_embedding(svd_list, common_1 = T, common_2 = F, distinct_1 = T, distinct_2 = F, 
                              vis_param = dcca_decomp$vis_param, noise_val = noise_val)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Distinct view",
     xlab = "UMAP 1", ylab = "UMAP 2")


# res <- custom_explained_variance(dcca_decomp)
# weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
# par(mfrow = c(1,1), mar = c(4, 4, 3.5, 0.5))
# information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F,
#                  main = "Simulation 2 (L2 Ratio)")
