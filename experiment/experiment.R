rm(list=ls())
set.seed(5)
n <- 100; K <- 2
common_space1 <- MASS::mvrnorm(n = n/2, mu = rep(0,K), Sigma = diag(K))
common_space2 <- MASS::mvrnorm(n = n/2, mu = rep(0,K), Sigma = diag(K))+20
common_space <- rbind(common_space1, common_space2)

p1 <- 5; p2 <- 10
transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)

mat_1 <- common_space %*% transform_mat_1 + scale(MASS::mvrnorm(n = n, mu = rep(0,p1), Sigma = diag(p1)), center = T, scale = F)
mat_2 <- common_space %*% transform_mat_2 + scale(MASS::mvrnorm(n = n, mu = rep(0,p2), Sigma = diag(p2)), center = T, scale = F)

dcca_obj <- dcca_factor(mat_1, mat_2, dims_1 = 1:K, dims_2 = 1:K, verbose = F)

membership_vec1 <- as.factor(rep(c("a","b"), each = n/2))
list_g1 <- construct_frnn(dcca_obj, data_1 = T, data_2 = F, nn = 25, membership_vec = membership_vec1,
                         verbose = F, bool_matrix = T)
res1 <- clisi_information(list_g1$c_g, list_g1$d_g, list_g1$e_g, membership_vec1,
                          verbose = F)
list_g2 <- construct_frnn(dcca_obj, data_1 = F, data_2 = T, nn = 25, membership_vec = membership_vec1,
                          verbose = F, bool_matrix = T)
res2 <- clisi_information(list_g2$c_g, list_g2$d_g, list_g2$e_g, membership_vec1,
                          verbose = F)

zz <- plot_clisi(res1, res2)
tmp2 <- cowplot::plot_grid(zz[[1]], zz[[2]])
tmp2
