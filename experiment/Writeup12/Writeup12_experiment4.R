rm(list=ls())

# high distinct 1, high distinct 2 with 3 clusters
set.seed(10)
n_clust <- 100
high <- 0.9; low <- 0.05
B_mat1 <- matrix(c(high, high, low, low, 
                   high, high, low, low, 
                   low, low, high, low,
                   low, low, low, high), 4, 4, byrow = T)
B_mat2 <- matrix(c(high, low, low, low,
                   low, high, high, low,
                   low, high, high, low,
                   low, low, low, high), 4, 4, byrow = T)
K <- ncol(B_mat1)
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust), rep(4, n_clust))
n <- length(membership_vec); true_membership_vec <- membership_vec
svd_u_1 <- generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)[,1:2]
svd_u_2 <- generate_sbm_orthogonal(B_mat2, membership_vec, centered = T)[,1:2]

set.seed(10)
p_1 <- 20; p_2 <- 40
svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
svd_v_1 <- generate_random_orthogonal(p_1, 2)
svd_v_2 <- generate_random_orthogonal(p_2, 2)

set.seed(10)
dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_1)

set.seed(10)
K <- 2
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec)
plot_scores(dcca_decomp, membership_vec = true_membership_vec, decomposition = T)
