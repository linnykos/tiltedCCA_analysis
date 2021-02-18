rm(list=ls())
library(Seurat)
source("../multiomicCCA_analysis/experiment/Writeup11/Writeup11_simulation_functions.R")

# let's start with a simple SBM, mostly common-space embedding
set.seed(10)
n_clust <- 100
B_mat <- matrix(c(0.9, 0.4, 0.1, 
                0.4, 0.9, 0.1,
                0.1, 0.1, 0.5), 3, 3)
K <- ncol(B_mat)
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
n <- length(membership_vec)
rho <- 1
prob_mat <- .compute_prob_mat(B_mat, membership_vec)
adj_mat <- .generate_adjaceny_mat(prob_mat)
tmp <- svd(crossprod(adj_mat))
score_1 <- adj_mat %*% tmp$u[,1:K]

set.seed(10)
p_1 <- 20; p_2 <- 40
coef_mat_1 <- matrix(stats::rnorm(K*p_1), K, p_1)
coef_mat_2 <- matrix(stats::rnorm(K*p_2), K, p_2)
coef_mat_1 <- coef_mat_1/max(svd(coef_mat_1)$d[1])*0.9
coef_mat_2 <- coef_mat_2/max(svd(coef_mat_2)$d[1])*0.9

set.seed(10)
dat <- generate_data_dmca(score_1, score_1, coef_mat_1, coef_mat_2)
#png("../../out/simulation/Writeup11/Writeup11_simulation1_true_score.png", height = 1200, width = 2500, units = "px", res = 300)
plot_scores(dat)
#graphics.off()

# png("../../out/simulation/Writeup11/Writeup11_simulation1_true_embedding.png", height = 800, width = 2000, units = "px", res = 300)
plot_embeddings(dat)
# graphics.off()

# try DCCA

set.seed(10)
dmca_res <- dmca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
res <- dmca_decomposition(dmca_res, verbose = F)

#png("../../out/simulation/Writeup11/Writeup11_simulation1_dcca.png", height = 800, width = 2000, units = "px", res = 300)
plot_embeddings(res)
#graphics.off()

####################################################
####################################################

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


svd_1 <- .svd_truncated(prob_mat1, K = 3); svd_2 <- .svd_truncated(prob_mat2, K = 3)
zz <- crossprod(prob_mat1, prob_mat2); Matrix::rankMatrix(zz); image(t(zz))
image(.compute_cca_aggregate_matrix(svd_1, svd_2))
tmp <- svd(zz)
score_1 <- prob_mat1 %*% tmp$u[,1:K]
score_2 <- prob_mat2 %*% tmp$v[,1:K]

set.seed(10)
p_1 <- 20; p_2 <- 40
coef_mat_1 <- matrix(stats::rnorm(K*p_1), K, p_1)
coef_mat_2 <- matrix(stats::rnorm(K*p_2), K, p_2)
coef_mat_1 <- coef_mat_1/max(svd(coef_mat_1)$d[1])*0.9
coef_mat_2 <- coef_mat_2/max(svd(coef_mat_2)$d[1])*0.9

set.seed(10)
dat <- generate_data_dmca(score_1, score_2, coef_mat_1, coef_mat_2)
plot_scores(dat)

plot_embeddings(dat, membership_vec = true_membership_vec)

set.seed(10)
dmca_res <- dmca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
res <- dmca_decomposition(dmca_res, verbose = F)

#png("../../out/simulation/Writeup11/Writeup11_simulation1_dcca.png", height = 800, width = 2000, units = "px", res = 300)
plot_embeddings(res, membership_vec = true_membership_vec)


####################################################
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
dat <- generate_data(score_1, score_2, coef_mat_1, coef_mat_2, reorthogonalize = T)

zlim <- range(c(dat$common_score, dat$distinct_score_1, dat$distinct_score_2))
c1 <- svd(dat$common_score)$d[1]; d1 <- svd(dat$distinct_score_1)$d[1]; d2 <- svd(dat$distinct_score_2)$d[1]
png("../../out/simulation/Writeup11/Writeup11_simulation3_dmca_true_score.png", height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,3))
image(t(dat$common_score), zlim = zlim, main = "Common score")
image(t(dat$distinct_score_1), zlim = zlim, main = paste0("Distinct score 1 (", round(d1/(c1+d1),2), "% var.)"))
image(t(dat$distinct_score_2), zlim = zlim, main = paste0("Distinct score 2 (", round(d2/(c1+d2),2), "% var.)"))
graphics.off()

# generate true images
true_dat <- construct_noiseless_data(dat$common_score, dat$distinct_score_1, dat$distinct_score_2,
                                     coef_mat_1, coef_mat_2)

png("../../out/simulation/Writeup11/Writeup11_simulation3_dmca_true_embedding.png", height = 800, width = 2000, units = "px", res = 300)
par(mfrow = c(1,3))
set.seed(10)
tmp <- extract_embedding(true_dat, common_1 = T, common_2 = T, distinct_1 = F, distinct_2 = F)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Common view (Truth)",
     xlab = "UMAP 1", ylab = "UMAP 2")
set.seed(10)
tmp <- extract_embedding(true_dat, common_1 = F, common_2 = F, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Distinct views (Truth)",
     xlab = "UMAP 1", ylab = "UMAP 2")
set.seed(10)
tmp <- extract_embedding(true_dat, common_1 = T, common_2 = T, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Entire view (Truth)",
     xlab = "UMAP 1", ylab = "UMAP 2")
graphics.off()

# try DCCA

set.seed(10)
dmca_res <- dmca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
res <- dmca_decomposition(dmca_res, verbose = F)

png("../../out/simulation/Writeup11/Writeup11_simulation3_dmca.png", height = 800, width = 2000, units = "px", res = 300)
par(mfrow = c(1,3))
set.seed(10)
tmp <- extract_embedding(res, common_1 = T, common_2 = T, distinct_1 = F, distinct_2 = F)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Common view (DCCA)",
     xlab = "UMAP 1", ylab = "UMAP 2")
set.seed(10)
tmp <- extract_embedding(res, common_1 = F, common_2 = F, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Distinct views (DCCA)",
     xlab = "UMAP 1", ylab = "UMAP 2")
set.seed(10)
tmp <- extract_embedding(res, common_1 = T, common_2 = T, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Entire view (DCCA)",
     xlab = "UMAP 1", ylab = "UMAP 2")
graphics.off()

