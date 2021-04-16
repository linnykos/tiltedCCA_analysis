# easy setting: low distinct 1, low distinct 2 with 3 clusters
rm(list=ls())
library(Seurat)
source("../multiomicCCA_analysis/experiment/Writeup12/Writeup12_simulation_functions.R")

set.seed(10)
n_clust <- 100
B_mat <- matrix(c(0.9, 0.2, 0.1, 
                  0.2, 0.9, 0.1,
                  0.1, 0.1, 0.5), 3, 3, byrow = T)
K <- ncol(B_mat); rho <- 1
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
n <- length(membership_vec); true_membership_vec <- membership_vec
svd_u_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)
svd_u_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)

set.seed(10)
p_1 <- 20; p_2 <- 40
svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
svd_v_1 <- generate_random_orthogonal(p_1, K-1)
svd_v_2 <- generate_random_orthogonal(p_2, K-1)

set.seed(10)
dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2)

set.seed(10)
K <- ncol(svd_u_1)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
dcca_decomp$common_perc

png("../../out/simulation/Writeup12/Writeup12_simulation1_data_umap.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation1_data_pca.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec, pca = T)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation1_scoreheatmap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation1_score.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores(dcca_decomp, membership_vec = true_membership_vec, decomposition = T)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation1_decomp.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5), mfrow = c(1,2))
plot_decomposition_2d(dcca_decomp$score_1[,1], dcca_decomp$score_2[,1], dcca_decomp$common_score[,1], 
                      main = bquote(atop(bold("Canonical dim. 1"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[1],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$common_perc[1],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")
plot_decomposition_2d(dcca_decomp$score_1[,2], dcca_decomp$score_2[,2], dcca_decomp$common_score[,2], 
                      main = bquote(atop(bold("Canonical dim. 2"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[2],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$common_perc[2],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation1_embedding1_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F, main_addition = "\nDataset 1")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation1_embedding2_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, main_addition = "\nDataset 2")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation1_embeddingboth_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T, main_addition = "\nBoth datasets")
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation1_embeddingboth_umap_nonoise.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T, add_noise = F, main_addition = "\nBoth datasets")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation1_embedding2_pca.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, pca = T, main_addition = "\nDataset 2")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation1_embedding1_pca.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F, pca = T, main_addition = "\nDataset 1")
graphics.off()
#######################################

# high distinct 1, low distinct 2 with 3 clusters
set.seed(10)
n_clust <- 100
B_mat1 <- matrix(c(0.9, 0, 0, 
                  0, 0.9, 0,
                  0, 0, 0.9), 3, 3, byrow = T)
B_mat2 <- matrix(c(0.9, 0.85, 0, 
                   0.85, 0.9, 0,
                   0, 0, 1), 3, 3, byrow = T)
K <- ncol(B_mat1)
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
n <- length(membership_vec); true_membership_vec <- membership_vec
svd_u_1 <- generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)
svd_u_2 <- generate_sbm_orthogonal(B_mat2, membership_vec, centered = T)

set.seed(10)
p_1 <- 20; p_2 <- 40
svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
svd_v_1 <- generate_random_orthogonal(p_1, K-1)
svd_v_2 <- generate_random_orthogonal(p_2, K-1)

set.seed(10)
dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2, noise_val = 0.1)
dat$common_perc

set.seed(10)
K <- ncol(svd_u_1)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
dcca_decomp$common_perc

# plot_embeddings(dcca_decomp, membership_vec = true_membership_vec)
png("../../out/simulation/Writeup12/Writeup12_simulation2_data_umap.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation2_data_pca.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec, pca = T)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation2_scoreheatmap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation2_score.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores(dcca_decomp, membership_vec = true_membership_vec, decomposition = T)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation2_decomp.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5), mfrow = c(1,2))
plot_decomposition_2d(dcca_decomp$score_1[,1], dcca_decomp$score_2[,1], dcca_decomp$common_score[,1], 
                      main = bquote(atop(bold("Canonical dim. 1"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[1],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$common_perc[1],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")
plot_decomposition_2d(dcca_decomp$score_1[,2], dcca_decomp$score_2[,2], dcca_decomp$common_score[,2], 
                      main = bquote(atop(bold("Canonical dim. 2"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[2],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$common_perc[2],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")
graphics.off()


png("../../out/simulation/Writeup12/Writeup12_simulation2_embedding1_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F, main_addition = "\nDataset 1")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation2_embedding2_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, main_addition = "\nDataset 2")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation2_embeddingboth_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T, main_addition = "\nBoth datasets")
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation2_embedding1_umap_nonoise.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F, main_addition = "\nDataset 1, No noise", add_noise = F)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation2_embedding2_umap_nonoise.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, main_addition = "\nDataset 2, No noise", add_noise = F)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation2_embeddingboth_umap_nonoise.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T, main_addition = "\nBoth datasets, No noise", add_noise = F)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation2_embedding2_umap_nonoise.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, add_noise = F, main_addition = "\nDataset 2")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation2_embedding2_pca.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, pca = T, main_addition = "\nDataset 2")
graphics.off()

##################################

# even more dramatic high distinct 1, low distinct 2 with 3 clusters
set.seed(10)
n_clust <- 100
high <- 0.9; low <- 0.05
B_mat1 <- matrix(c(0.9, 0.1, 0.1,
                   0.1, 0.9, 0.1,
                   0.1, 0.1, 0.9), 3, 3, byrow = T)
K <- ncol(B_mat1)
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
n <- length(membership_vec); true_membership_vec <- membership_vec
svd_u_1 <- generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)[,1:2]
svd_u_2 <- generate_random_orthogonal(n, 2, centered = T)

set.seed(10)
p_1 <- 20; p_2 <- 40
svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
svd_v_1 <- generate_random_orthogonal(p_1, 2)
svd_v_2 <- generate_random_orthogonal(p_2, 2)

set.seed(10)
dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_1, noise_val = 2)

set.seed(10)
K <- 2
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
dcca_decomp$common_perc

png("../../out/simulation/Writeup12/Writeup12_simulation3_data_umap.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation3_data_pca.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec, pca = T)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation3_score_data.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()


png("../../out/simulation/Writeup12/Writeup12_simulation3_scoreheatmap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation3_score.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores(dcca_decomp, membership_vec = true_membership_vec, decomposition = T)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation3_decomp.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5), mfrow = c(1,2))
plot_decomposition_2d(dcca_decomp$score_1[,1], dcca_decomp$score_2[,1], dcca_decomp$common_score[,1], 
                      main = bquote(atop(bold("Canonical dim. 1"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[1],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$common_perc[1],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")
plot_decomposition_2d(dcca_decomp$score_1[,2], dcca_decomp$score_2[,2], dcca_decomp$common_score[,2], 
                      main = bquote(atop(bold("Canonical dim. 2"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[2],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$common_perc[2],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")
graphics.off()


png("../../out/simulation/Writeup12/Writeup12_simulation3_embedding1_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F, main_addition = "\nDataset 1")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation3_embedding2_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, main_addition = "\nDataset 2")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation3_embeddingboth_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T, main_addition = "\nBoth datasets")
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation3_embedding1_umap_nonoise.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F, main_addition = "\nDataset 1, No noise", add_noise = F)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation3_embedding2_umap_nonoise.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, main_addition = "\nDataset 2, No noise", add_noise = F)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation3_embeddingboth_umap_nonoise.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T, main_addition = "\nBoth datasets, No noise", add_noise = F)
graphics.off()

set.seed(10)
K <- ncol(svd_u_1)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F, fix_common_perc = T)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
png("../../out/simulation/Writeup12/Writeup12_simulation3_scoreheatmap_fixedcommon.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation3_score_fixedcommon.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores(dcca_decomp, membership_vec = true_membership_vec, decomposition = T)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation3_decomp_fixedcommon.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5), mfrow = c(1,2))
plot_decomposition_2d(dcca_decomp$score_1[,1], dcca_decomp$score_2[,1], dcca_decomp$common_score[,1], 
                      main = "Canonical dim. 1", xlab = "Dataset 1", ylab = "Dataset 2")
plot_decomposition_2d(dcca_decomp$score_1[,2], dcca_decomp$score_2[,2], dcca_decomp$common_score[,2], 
                      main = "Canonical dim. 2", xlab = "Dataset 1", ylab = "Dataset 2")
graphics.off()

plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec)
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F, pca = T)
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, pca = T)
plot_scores(dcca_decomp, membership_vec = true_membership_vec, decomposition = T)
par(mfrow = c(1,2))
plot_decomposition_2d(dcca_decomp$score_1[,1], dcca_decomp$score_2[,1], dcca_decomp$common_score[,1])
plot_decomposition_2d(dcca_decomp$score_1[,2], dcca_decomp$score_2[,2], dcca_decomp$common_score[,2])

####################################3

# high distinct 1, high distinct 2 with 4 clusters. 
# distinct 1 separates (1) from (2) and 
# distinct 2 seperates (3) from (4) analogously

set.seed(10)
B_mat <- matrix(c(0.9, 0.4, 0.1,
                  0.4, 0.9, 0.1,
                  0.1, 0.1, 0.3), 3, 3, byrow = T)
K <- ncol(B_mat); n_clust <- 100

true_membership_vec <- rep(1:4, each = n_clust)
n <- length(true_membership_vec)

membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, 2*n_clust))
svd_u_1 <- generate_sbm_orthogonal(B_mat, membership_vec, centered = T)
membership_vec <- c(rep(3, 2*n_clust), rep(1, n_clust), rep(2, n_clust))
svd_u_2 <- generate_sbm_orthogonal(B_mat, membership_vec, centered = T)

set.seed(10)
p_1 <- 20; p_2 <- 40
svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
svd_v_1 <- generate_random_orthogonal(p_1, K-1)
svd_v_2 <- generate_random_orthogonal(p_2, K-1)

set.seed(10)
dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2)

set.seed(10)
K <- ncol(svd_u_1)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
dcca_decomp$common_perc

png("../../out/simulation/Writeup12/Writeup12_simulation4_data_umap.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation4_data_pca.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec, pca = T)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation4_scoreheatmap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation4_score.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores(dcca_decomp, membership_vec = true_membership_vec, decomposition = T)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation4_decomp.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5), mfrow = c(1,2))
plot_decomposition_2d(dcca_decomp$score_1[,1], dcca_decomp$score_2[,1], dcca_decomp$common_score[,1], 
                      main = bquote(atop(bold("Canonical dim. 1"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[1],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$common_perc[1],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")
plot_decomposition_2d(dcca_decomp$score_1[,2], dcca_decomp$score_2[,2], dcca_decomp$common_score[,2], 
                      main = bquote(atop(bold("Canonical dim. 2"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[2],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$common_perc[2],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")
graphics.off()


png("../../out/simulation/Writeup12/Writeup12_simulation4_embedding1_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F, main_addition = "\nDataset 1")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation4_embedding2_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, main_addition = "\nDataset 2")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation4_embeddingboth_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T, main_addition = "\nBoth datasets")
graphics.off()


plot_scores(dcca_decomp, true_membership_vec)
plot_scores(dcca_decomp, true_membership_vec, decomposition = T)

##################################

# high distinct 1, high distinct 2 with 5 clusters
# the common loading separates (1234) from (5)

set.seed(10)
B_mat <- matrix(c(0.9, 0.4, 0.1, 
                  0.4, 0.9, 0.1, 
                  0.1, 0.1, 0.9), 3, 3, byrow = T)
K <- ncol(B_mat); n_clust <- 100
true_membership_vec <- rep(1:5, each = n_clust)
n <- length(true_membership_vec)

membership_vec <-  c(rep(1, 2*n_clust), rep(2, 2*n_clust), rep(3, n_clust))
svd_u_1 <- generate_sbm_orthogonal(B_mat, membership_vec, centered = T)
membership_vec <-  c(rep(1, n_clust), rep(2, n_clust), rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
svd_u_2 <- generate_sbm_orthogonal(B_mat, membership_vec, centered = T)

set.seed(10)
p_1 <- 40; p_2 <- 40
svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
svd_v_1 <- generate_random_orthogonal(p_1, K-1)
svd_v_2 <- generate_random_orthogonal(p_2, K-1)

set.seed(10)
dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2,
                     noise_val = 2)
dat$common_perc

set.seed(10)
K <- ncol(svd_u_1)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F, num_neigh = 80)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
dcca_decomp$common_perc

png("../../out/simulation/Writeup12/Writeup12_simulation5_data_umap.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation5_data_pca.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec, pca = T)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation5_scoreheatmap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation5_score.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores(dcca_decomp, membership_vec = true_membership_vec, decomposition = T)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation5_decomp.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5), mfrow = c(1,2))
plot_decomposition_2d(dcca_decomp$score_1[,1], dcca_decomp$score_2[,1], dcca_decomp$common_score[,1], 
                      main = bquote(atop(bold("Canonical dim. 1"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[1],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$common_perc[1],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")
plot_decomposition_2d(dcca_decomp$score_1[,2], dcca_decomp$score_2[,2], dcca_decomp$common_score[,2], 
                      main = bquote(atop(bold("Canonical dim. 2"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[2],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$common_perc[2],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation5_embedding1_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F, main_addition = "\nDataset 1")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation5_embedding2_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, main_addition = "\nDataset 2")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation5_embeddingboth_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T, main_addition = "\nBoth datasets")
graphics.off()

set.seed(10)
K <- 2
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F, fix_common_perc = T)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
dcca_decomp$common_perc
plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec)
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F, pca = T)
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, pca = T)
plot_scores(dcca_decomp, membership_vec = true_membership_vec, decomposition = T)

##################################

# same setting as above, but now w/ different number variables
set.seed(10)
B_mat <- matrix(c(0.9, 0.4, 0.1, 
                  0.4, 0.9, 0.1, 
                  0.1, 0.1, 0.9), 3, 3, byrow = T)
K <- ncol(B_mat); n_clust <- 100
true_membership_vec <- rep(1:5, each = n_clust)
n <- length(true_membership_vec)

membership_vec <-  c(rep(1, 2*n_clust), rep(2, 2*n_clust), rep(3, n_clust))
svd_u_1 <- generate_sbm_orthogonal(B_mat, membership_vec, centered = T)
membership_vec <-  c(rep(1, n_clust), rep(2, n_clust), rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
svd_u_2 <- generate_sbm_orthogonal(B_mat, membership_vec, centered = T)

set.seed(10)
p_1 <- 20; p_2 <- 40
svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_1)*c(0.75,0.5)
# TAKEAWAY: It might be fine, but both matrices need to have the same 'amount of signal'
#  via the singular values...
svd_v_1 <- generate_random_orthogonal(p_1, K-1)
svd_v_2 <- generate_random_orthogonal(p_2, K-1)

set.seed(10)
dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2,
                     noise_val = 2)
dat$common_perc

set.seed(10)
K <- ncol(svd_u_1)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F, num_neigh = 80)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)
dcca_decomp$common_perc

png("../../out/simulation/Writeup12/Writeup12_simulation6_data_umap.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation6_data_pca.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec, pca = T)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation6_scoreheatmap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation6_score.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_scores(dcca_decomp, membership_vec = true_membership_vec, decomposition = T)
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation6_decomp.png", height = 1000, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5), mfrow = c(1,2))
plot_decomposition_2d(dcca_decomp$score_1[,1], dcca_decomp$score_2[,1], dcca_decomp$common_score[,1], 
                      main = bquote(atop(bold("Canonical dim. 1"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[1],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$common_perc[1],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")
plot_decomposition_2d(dcca_decomp$score_1[,2], dcca_decomp$score_2[,2], dcca_decomp$common_score[,2], 
                      main = bquote(atop(bold("Canonical dim. 2"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[2],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$common_perc[2],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation6_embedding1_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F, main_addition = "\nDataset 1")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation6_embedding2_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, main_addition = "\nDataset 2")
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation6_embeddingboth_umap.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T, main_addition = "\nBoth datasets")
graphics.off()

png("../../out/simulation/Writeup12/Writeup12_simulation6_embedding1_umap_nonoise.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F, main_addition = "\nDataset 1, No noise", add_noise = F)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation6_embedding2_umap_nonoise.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, main_addition = "\nDataset 2, No noise", add_noise = F)
graphics.off()
png("../../out/simulation/Writeup12/Writeup12_simulation6_embeddingboth_umap_nonoise.png", height = 700, width = 1700, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T, main_addition = "\nBoth datasets, No noise", add_noise = F)
graphics.off()


plot_scores(dcca_decomp, membership_vec = true_membership_vec)
plot_scores(dcca_decomp, membership_vec = true_membership_vec, decomposition = T)
plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec)

#3333333###333333333################

a <- 2
b_seq <- seq(0.01, 5, length.out = 100)
y_seq <- sapply(b_seq, function(b){.sigmoid_ratio(a,b)})
plot(b_seq-a, y_seq, xlab = "(b-a) for a=2", ylab = "hinge(a,b)")
lines(rep(0,2), c(-10,10), col = 'red', lty = 2)
lines(c(-10,10), rep(0.5,2), col = 'red', lty = 2)
points(b_seq-a, y_seq, pch = 16)
