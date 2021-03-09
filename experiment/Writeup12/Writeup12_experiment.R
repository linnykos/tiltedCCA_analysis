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
                          noise_func = function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat, sd = 0.1), nrow(mat), ncol(mat))})
plot_data(dat, membership_vec = true_membership_vec, raw = T, K = 3)

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)


plot_embeddings(dcca_decomp, membership_vec = true_membership_vec)
plot_data(dcca_decomp, membership_vec = true_membership_vec)

#########

n <- nrow(dat$mat_1)
cs <- dcca_res$common_score
full_mat <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2
zz <- (diag(n) - cs %*% solve(crossprod(cs)) %*% t(cs)) %*% full_mat

quantile(full_mat)
quantile(zz)

# plot some PCAs
par(mfrow = c(1,2))
full_mat1 <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
pc1 <- prcomp(full_mat1)
plot(pc1$x[,1], pc1$x[,2], col = true_membership_vec, pch = 16, asp = T)
full_mat2 <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2
pc2 <- prcomp(full_mat2)
plot(pc2$x[,1], pc2$x[,2], col = true_membership_vec, pch = 16, asp = T)


pc_c<- prcomp(dcca_decomp$common_score)
plot(pc_c$x[,1], pc_c$x[,2], col = true_membership_vec, pch = 16, asp = T)

pc3 <- prcomp(dcca_decomp$common_mat_2)
plot(pc3$x[,1], pc3$x[,2], col = true_membership_vec, pch = 16, asp = T)


par(mfrow = c(1,2))
plot(NA, xlim = range( dcca_res$score_1), ylim = c(0.5,3.5))
n <- nrow(dcca_res$score_1)
for(i in 1:3){
  points(x = dcca_res$score_1[,i], y = runif(n, min = i-.2, max = i+.2), col = true_membership_vec, pch = 16)
}
plot(NA, xlim = range( dcca_res$score_2), ylim = c(0.5,3.5))
n <- nrow(dcca_res$score_2)
for(i in 1:3){
  points(x = dcca_res$score_2[,i], y = runif(n, min = i-.2, max = i+.2), col = true_membership_vec, pch = 16)
}

par(mfrow = c(1,1))
plot(dcca_res$score_1[,2], dcca_res$score_2[,2], pch = 16, col = rep(1:3, each = 100), asp = T)
i <- 1
lines(rep(dcca_res$score_1[i,2], 2), c(-10,10), col = "red")
lines( c(-10,10), rep(dcca_res$score_2[i,2], 2), col = "red")

par(mfrow = c(1,3))
for(i in 1:3){
  plot(dcca_res$score_1[,i], dcca_res$score_2[,i], pch = 16, col = rep(1:3, each = 100), asp = T)
  
}


###############
nn_1 <- RANN::nn2(dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1, k = 50)
nn_2 <- RANN::nn2(dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2, k = 50)

for(i in 1:3){
  print(.latent_common_perc(dcca_res$score_1[,i], dcca_res$score_2[,i], nn_1, nn_2))
}
