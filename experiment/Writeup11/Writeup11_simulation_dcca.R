rm(list=ls())
library(Seurat)
source("../multiomicCCA_analysis/experiment/Writeup11/Writeup11_simulation_functions.R")

set.seed(10)
n_clust <- 100
B_mat <- matrix(c(0.9, 0.4, 0.1, 
                  0.4, 0.9, 0.1,
                  0.1, 0.1, 0.5), 3, 3)
K <- ncol(B_mat)
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
n <- length(membership_vec)
true_membership_vec <- membership_vec
rho <- 1
score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)

set.seed(10)
p_1 <- 20; p_2 <- 40
coef_mat_1 <- matrix(stats::rnorm(K*p_1), K, p_1)
coef_mat_2 <- matrix(stats::rnorm(K*p_2), K, p_2)

set.seed(10)
dat <- generate_data_dcca(score_1, score_2, coef_mat_1, coef_mat_2)
plot_scores(dat, mode = 1)
plot_embeddings(dat, membership_vec = true_membership_vec, mode = "dcca")

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
plot_embeddings(dcca_decomp, membership_vec = true_membership_vec, mode = "dcca")

dcca_decomp2 <- dcca_variance_decomposition(dcca_res, rank_c = K, verbose = F)
res <- explained_variance(dcca_decomp2, verbose = F)
weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
par(mfrow = c(1,1), mar = c(6, 4, 3.5, 0.5))
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)

tmp <- dcca_information_weight(mat = dcca_decomp$common_mat_1+dcca_decomp$distinct_mat_1,
                               common_mat = dcca_decomp$common_mat_1)
weight_vec1 <- tmp$alpha_vec
tmp <- dcca_information_weight(mat = dcca_decomp$common_mat_2+dcca_decomp$distinct_mat_2,
                               common_mat = dcca_decomp$common_mat_2)
weight_vec2 <- tmp$alpha_vec
weight_mat <- data.frame(Mode_1 = 1-weight_vec1, Mode_2 = 1-weight_vec2, cell_type = true_membership_vec)
par(mfrow = c(1,1), mar = c(6, 4, 3.5, 0.5))
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)

#####################################

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
dat <- generate_data_dcca(score_1, score_2, coef_mat_1, coef_mat_2); class(dat) <- "dcca_decomp"
plot_scores(dat, mode = 1)
plot_embeddings(dat, membership_vec = true_membership_vec, mode = "dcca")
pipeline_information_plot(dat, true_membership_vec)
res <- custom_explained_variance2(dat, true_membership_vec)
weight_mat <- data.frame(Mode_1 = res$weight_vec_1, Mode_2 = res$weight_vec_2, cell_type = 1:length(unique(true_membership_vec)))
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), 
                 reorder_types = F, plot_individual = F)

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
plot_embeddings(dcca_decomp, membership_vec = true_membership_vec, mode = "dcca")
pipeline_information_plot(dcca_decomp, true_membership_vec)

dcca_decomp2 <- dcca_variance_decomposition(dcca_res, rank_c = K, verbose = F)
res <- explained_variance(dcca_decomp2, verbose = F)
weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
par(mfrow = c(1,1), mar = c(6, 4, 3.5, 0.5))
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)

# try quadratic
res <- explained_variance(dcca_decomp2, weight_func = function(x){x^2}, verbose = F)
weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)
# try exponential
res <- explained_variance(dcca_decomp2, weight_func = function(x){exp(x)}, verbose = F)
weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)
# try sum-of-square partition
res <- custom_explained_variance(dcca_decomp)
weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)
# try sum-of-square partition squared
res <- custom_explained_variance(dcca_decomp, weight_func = function(x){x^2})
weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)
# try sum-of-square partition exp
res <- custom_explained_variance(dcca_decomp, weight_func = function(x){exp(x)})
weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)


res <- custom_explained_variance2(dcca_decomp, true_membership_vec)
weight_mat <- data.frame(Mode_1 = res$weight_vec_1, Mode_2 = res$weight_vec_2, cell_type = 1:length(unique(true_membership_vec)))
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), 
                 reorder_types = F, plot_individual = F)

####################################################

rm(list=ls())
library(Seurat)
source("../multiomicCCA_analysis/experiment/Writeup11/Writeup11_simulation_functions.R")

# try a setting where the common loading separates (1234) from (5), 
# distinct 1 separates (12) from (34), and distinct 2 separates (13) from (24)
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
dat <- generate_data_dcca(score_1, score_2, coef_mat_1, coef_mat_2)
plot_scores(dat, mode = 1)
plot_embeddings(dat, membership_vec = true_membership_vec, mode = "dcca")

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
plot_embeddings(dcca_decomp, membership_vec = true_membership_vec, mode = "dcca")

tmp <- dcca_information_weight(mat = dcca_decomp$common_mat_1+dcca_decomp$distinct_mat_1,
                               common_mat = dcca_decomp$common_mat_1)
weight_vec1 <- tmp$alpha_vec
tmp <- dcca_information_weight(mat = dcca_decomp$common_mat_2+dcca_decomp$distinct_mat_2,
                               common_mat = dcca_decomp$common_mat_2)
weight_vec2 <- tmp$alpha_vec
weight_mat <- data.frame(Mode_1 = 1-weight_vec1, Mode_2 = 1-weight_vec2, cell_type = true_membership_vec)
par(mfrow = c(1,1), mar = c(6, 4, 3.5, 0.5))
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)

dcca_decomp2 <- dcca_variance_decomposition(dcca_res, rank_c = K, verbose = F)
res <- explained_variance(dcca_decomp2, verbose = F)
weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
par(mfrow = c(1,1), mar = c(6, 4, 3.5, 0.5))
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)

# try quadratic
res <- explained_variance(dcca_decomp2, weight_func = function(x){x^2}, verbose = F)
weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)
# try exponential
res <- explained_variance(dcca_decomp2, weight_func = function(x){exp(x)}, verbose = F)
weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)
# try sum-of-square partition
res <- custom_explained_variance(dcca_decomp)
weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)
# try sum-of-square partition squared
res <- custom_explained_variance(dcca_decomp, weight_func = function(x){x^2})
weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)
# try sum-of-square partition exp
res <- custom_explained_variance(dcca_decomp, weight_func = function(x){exp(x)})
weight_mat <- data.frame(Mode_1 = res$cell_weight_vec1, Mode_2 = res$cell_weight_vec2, cell_type = true_membership_vec)
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F)

res <- custom_explained_variance2(dcca_decomp, true_membership_vec)
weight_mat <- data.frame(Mode_1 = res$weight_vec_1, Mode_2 = res$weight_vec_2, cell_type = 1:length(unique(true_membership_vec)))
information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), 
                 reorder_types = F, plot_individual = F)


