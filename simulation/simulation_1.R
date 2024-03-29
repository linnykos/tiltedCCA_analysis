rm(list=ls())
library(tiltedCCA)
set.seed(10)

################################
# Step 1: Generate the data, where Modality 1 is informative, but not Modality 2
################################

n_clust <- 100
B_mat1 <- matrix(c(0.9, 0.1, 0.1,
                   0.1, 0.9, 0.1,
                   0.1, 0.1, 0.9), 3, 3, byrow = T)
K <- ncol(B_mat1)
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
true_cluster <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
n <- length(membership_vec); true_membership_vec <- membership_vec
svd_u_1 <- tiltedCCA::generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)[,1:2]
svd_u_2 <- tiltedCCA::generate_random_orthogonal(n, 2, centered = T)

p_1 <- 10; p_2 <- 10
svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
svd_v_1 <- tiltedCCA::generate_random_orthogonal(p_1, 2)
svd_v_2 <- tiltedCCA::generate_random_orthogonal(p_2, 2)

mat_1 <- tcrossprod(svd_u_1 %*% diag(svd_d_1), svd_v_1)
mat_2 <- tcrossprod(svd_u_2 %*% diag(svd_d_2), svd_v_2)

clustering_1 <- factor(stats::kmeans(mat_1, centers = 3)$cluster)
clustering_2 <- factor(rep(1, length(membership_vec)))

rownames(mat_1) <- paste0("n", 1:nrow(mat_1))
rownames(mat_2) <- paste0("n", 1:nrow(mat_2))
colnames(mat_1) <- paste0("g", 1:ncol(mat_1))
colnames(mat_2) <- paste0("p", 1:ncol(mat_2))

################################
# Step 2: Apply Tilted-CCA
################################

set.seed(10)
multiSVD_obj <- tiltedCCA::create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                            dims_1 = 1:2, dims_2 = 1:2,
                                            center_1 = F, center_2 = F,
                                            normalize_row = T,
                                            normalize_singular_value = F,
                                            recenter_1 = F, recenter_2 = F,
                                            rescale_1 = F, rescale_2 = F,
                                            scale_1 = F, scale_2 = F)
multiSVD_obj <- tiltedCCA::form_metacells(input_obj = multiSVD_obj,
                                           large_clustering_1 = clustering_1, 
                                           large_clustering_2 = clustering_2,
                                           num_metacells = NULL)
multiSVD_obj <- tiltedCCA::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 2,
                                         num_neigh = 10,
                                         bool_cosine = T,
                                         bool_intersect = T,
                                         min_deg = 1)
multiSVD_obj <- tiltedCCA::tiltedCCA(input_obj = multiSVD_obj)
multiSVD_obj <- tiltedCCA::fine_tuning(input_obj = multiSVD_obj,
                                        verbose = 0)
multiSVD_obj <- tiltedCCA::tiltedCCA_decomposition(multiSVD_obj)

################################
# Step 3: Plot the data
################################

# png("simulation1_data.png", height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(multiSVD_obj$svd_1$u[,1], multiSVD_obj$svd_1$u[,2],
     main = "Modality 1",
     xlab = "PCA's dim. 1", ylab = "PCA's dim. 2",
     pch = 16, col = true_cluster, asp = T)
plot(multiSVD_obj$svd_2$u[,1], multiSVD_obj$svd_2$u[,2],
     main = "Modality 2",
     xlab = "PCA's dim. 1", ylab = "PCA's dim. 2",
     pch = 16, col = true_cluster, asp = T)
# graphics.off()

################################
# Step 4: Plot Tilted-CCA's result
################################

names(multiSVD_obj)

# png("simulation1_tcca.png", height = 1200, width = 3000, res = 300, units = "px")
par(mfrow = c(1,3))
plot(multiSVD_obj$tcca_obj$common_score[,1], multiSVD_obj$tcca_obj$common_score[,2],
     main = "Common embedding",
     xlab = "Common's dim. 1", ylab = "Common's dim. 2",
     pch = 16, col = true_cluster, asp = T)
plot(multiSVD_obj$tcca_obj$distinct_score_1[,1], multiSVD_obj$tcca_obj$distinct_score_1[,2],
     main = "Modality 1's distinct embedding",
     xlab = "Distinct-1's dim. 1", ylab = "Distinct-1's dim. 2",
     pch = 16, col = true_cluster, asp = T)
plot(multiSVD_obj$tcca_obj$distinct_score_2[,1], multiSVD_obj$tcca_obj$distinct_score_2[,2],
     main = "Modality 2's distinct embedding",
     xlab = "Distinct-2's dim. 1", ylab = "Distinct-2's dim. 2",
     pch = 16, col = true_cluster, asp = T)
# graphics.off()

################################
# Step 5: For comparison, plot Consensus PCA
################################

set.seed(10)
consensus_pca <- tiltedCCA:::consensus_pca(mat_1 = mat_1, mat_2 = mat_2,
                                           dims_1 = 1:2, dims_2 = 1:2,
                                           dims_consensus = 1:2,
                                           apply_pca = T,
                                           center_1 = F, center_2 = F,
                                           center_consensus = F,
                                           recenter_1 = F, recenter_2 = F,
                                           rescale_1 = F, rescale_2 = F,
                                           scale_1 = F, scale_2 = F,
                                           scale_consensus = F,
                                           verbose = 0)

# png("simulation1_consensuspca.png", height = 1200, width = 1200, res = 300, units = "px")
par(mfrow = c(1,1))
plot(consensus_pca$dimred_consensus[,1], consensus_pca$dimred_consensus[,2],
     main = "Consensus PCA embedding",
     xlab = "Consensus PCA's dim. 1", ylab = "Consensus PCA's dim. 2",
     pch = 16, col = true_cluster, asp = T)
# graphics.off()
