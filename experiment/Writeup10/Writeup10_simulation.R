rm(list=ls())
library(Seurat)
source("../multiomicCCA_analysis/experiment/Writeup10/Writeup10_simulation_functions.R")

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
score_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)
score_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec)

set.seed(10)
p_1 <- 20; p_2 <- 40
coef_mat_1 <- matrix(stats::rnorm(K*p_1), K, p_1)
coef_mat_2 <- matrix(stats::rnorm(K*p_2), K, p_2)

set.seed(10)
dat <- generate_data(score_1, score_2, coef_mat_1, coef_mat_2)

par(mfrow = c(1,3))
zlim <- range(c(dat$common_score, dat$distinct_score_1, dat$distinct_score_2))
image(t(dat$common_score), zlim = zlim)
image(t(dat$distinct_score_1), zlim = zlim)
image(t(dat$distinct_score_2), zlim = zlim)

# generate true images
true_dat <- construct_noiseless_data(dat$common_score, dat$distinct_score_1, dat$distinct_score_2,
                                     coef_mat_1, coef_mat_2)
par(mfrow = c(1,3))
set.seed(10)
tmp <- extract_embedding(true_dat, common = T, distinct_1 = F, distinct_2 = F)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Common view (Truth)")
set.seed(10)
tmp <- extract_embedding(true_dat, common = F, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Distinct view (Truth)")
set.seed(10)
tmp <- extract_embedding(true_dat, common = T, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Entire view (Truth)")


# try Seurat 

set.seed(10)
seurat_obj <- analyze_seurat_pipeline(dat$mat_1, dat$mat_2)

par(mfrow = c(1,3))
tmp <- seurat_obj[["umap1"]]@cell.embeddings
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Mode 1")
tmp <- seurat_obj[["umap2"]]@cell.embeddings
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Mode 2")
tmp <- seurat_obj[["wnn.umap"]]@cell.embeddings
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "WNN")

# try DCCA

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

par(mfrow = c(1,3))
set.seed(10)
tmp <- extract_embedding(res, common = T, distinct_1 = F, distinct_2 = F)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Common view (DCCA)")
set.seed(10)
tmp <- extract_embedding(res, common = F, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Distinct view (DCCA)")
set.seed(10)
tmp <- extract_embedding(res, common = T, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Entire view (DCCA)")


# concatenated PCA
svd_res_1 <- .svd_truncated(dat$mat_1, K)
svd_res_2 <- .svd_truncated(dat$mat_2, K)
zz2 <- cbind(.mult_mat_vec(svd_res_1$u, svd_res_1$d), 
            .mult_mat_vec(svd_res_2$u, svd_res_2$d))
set.seed(10)
tmp2 <- Seurat::RunUMAP(zz2, verbose = F)@cell.embeddings
par(mfrow = c(1,1))
plot(tmp2[,1], tmp2[,2], asp = T, pch = 16, col = membership_vec, main = "PCA concatenated")

####################################################
####################################################

rm(list=ls())
library(Seurat)
source("../multiomicCCA_analysis/experiment/Writeup10/Writeup10_simulation_functions.R")

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
dat <- generate_data(score_1, score_2, coef_mat_1, coef_mat_2)

par(mfrow = c(1,3))
zlim <- range(c(dat$common_score, dat$distinct_score_1, dat$distinct_score_2))
image(t(dat$common_score), zlim = zlim)
image(t(dat$distinct_score_1), zlim = zlim)
image(t(dat$distinct_score_2), zlim = zlim)

# generate true images
true_dat <- construct_noiseless_data(dat$common_score, dat$distinct_score_1, dat$distinct_score_2,
                                     coef_mat_1, coef_mat_2)
par(mfrow = c(1,3))
set.seed(10)
tmp <- extract_embedding(true_dat, common = T, distinct_1 = F, distinct_2 = F)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Common view (Truth)")
set.seed(10)
tmp <- extract_embedding(true_dat, common = F, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Distinct view (Truth)")
set.seed(10)
tmp <- extract_embedding(true_dat, common = T, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Entire view (Truth)")

# try Seurat 

set.seed(10)
seurat_obj <- analyze_seurat_pipeline(dat$mat_1, dat$mat_2)

par(mfrow = c(1,3))
tmp <- seurat_obj[["umap1"]]@cell.embeddings
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Mode 1")
tmp <- seurat_obj[["umap2"]]@cell.embeddings
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Mode 2")
tmp <- seurat_obj[["wnn.umap"]]@cell.embeddings
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "WNN")

# try DCCA

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

par(mfrow = c(1,3))
set.seed(10)
tmp <- extract_embedding(res, common = T, distinct_1 = F, distinct_2 = F)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Common view (DCCA)")
set.seed(10)
tmp <- extract_embedding(res, common = F, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Distinct view (DCCA)")
set.seed(10)
tmp <- extract_embedding(res)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Entire view (DCCA)")

# concatenated PCA
svd_res_1 <- .svd_truncated(dat$mat_1, K)
svd_res_2 <- .svd_truncated(dat$mat_2, K)
zz2 <- cbind(.mult_mat_vec(svd_res_1$u, svd_res_1$d), 
             .mult_mat_vec(svd_res_2$u, svd_res_2$d))
set.seed(10)
tmp2 <- Seurat::RunUMAP(zz2, verbose = F)@cell.embeddings
par(mfrow = c(1,1))
plot(tmp2[,1], tmp2[,2], asp = T, pch = 16, col = true_membership_vec, main = "PCA concatenated")


####################################################
####################################################

rm(list=ls())
library(Seurat)
source("../multiomicCCA_analysis/experiment/Writeup10/Writeup10_simulation_functions.R")

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
dat <- generate_data(score_1, score_2, coef_mat_1, coef_mat_2)

par(mfrow = c(1,3))
zlim <- range(c(dat$common_score, dat$distinct_score_1, dat$distinct_score_2))
image(t(dat$common_score), zlim = zlim)
image(t(dat$distinct_score_1), zlim = zlim)
image(t(dat$distinct_score_2), zlim = zlim)

# generate true images
true_dat <- construct_noiseless_data(dat$common_score, dat$distinct_score_1, dat$distinct_score_2,
                                     coef_mat_1, coef_mat_2)
par(mfrow = c(1,3))
set.seed(10)
tmp <- extract_embedding(true_dat, common = T, distinct_1 = F, distinct_2 = F)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Common view (Truth)")
set.seed(10)
tmp <- extract_embedding(true_dat, common = F, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Distinct view (Truth)")
set.seed(10)
tmp <- extract_embedding(true_dat, common = T, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Entire view (Truth)")


# set.seed(10)
# tmp <- Seurat::RunUMAP(common_loading, verbose = F)@cell.embeddings
# par(mfrow = c(1,1))
# plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "All DCCA")

# try Seurat 

set.seed(10)
seurat_obj <- analyze_seurat_pipeline(dat$mat_1, dat$mat_2)

par(mfrow = c(1,3))
tmp <- seurat_obj[["umap1"]]@cell.embeddings
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Mode 1")
tmp <- seurat_obj[["umap2"]]@cell.embeddings
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Mode 2")
tmp <- seurat_obj[["wnn.umap"]]@cell.embeddings
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "WNN")

# try DCCA

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
res <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

par(mfrow = c(1,3))
set.seed(10)
tmp <- extract_embedding(res, common = T, distinct_1 = F, distinct_2 = F)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Common view (DCCA)")
set.seed(10)
tmp <- extract_embedding(res, common = F, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Distinct view (DCCA)")
set.seed(10)
tmp <- extract_embedding(res, common = T, distinct_1 = T, distinct_2 = T)
plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = true_membership_vec, main = "Entire view (DCCA)")

# concatenated PCA
svd_res_1 <- .svd_truncated(dat$mat_1, K)
svd_res_2 <- .svd_truncated(dat$mat_2, K)
zz2 <- cbind(.mult_mat_vec(svd_res_1$u, svd_res_1$d), 
             .mult_mat_vec(svd_res_2$u, svd_res_2$d))
set.seed(10)
tmp2 <- Seurat::RunUMAP(zz2, verbose = F)@cell.embeddings
par(mfrow = c(1,1))
plot(tmp2[,1], tmp2[,2], asp = T, pch = 16, col = true_membership_vec, main = "PCA concatenated")


