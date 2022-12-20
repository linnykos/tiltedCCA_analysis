rm(list=ls())
library(tiltedCCA)
set.seed(10)

################################
# Step 1: Generate the data, where Modality 1 is informative, but not Modality 2
################################

n_each <- 100
true_cluster <- rep(1:3, each = n_each)
mat_1 <- do.call(rbind, lapply(1:3, function(i){
  if(i == 1){
    MASS::mvrnorm(n = n_each, mu = c(0,0,0), Sigma = diag(3))  
  } else if(i == 2) {
    MASS::mvrnorm(n = n_each, mu = c(15,0,0), Sigma = diag(3)) 
  } else {
    MASS::mvrnorm(n = n_each, mu = c(0,15,0), Sigma = diag(3)) 
  }
}))

mat_2 <- do.call(rbind, lapply(1:3, function(i){
  if(i %in% c(1,3)){
    MASS::mvrnorm(n = n_each, mu = c(0,0,0), Sigma = diag(3)) 
  } else {
    MASS::mvrnorm(n = n_each, mu = c(6,6,0), Sigma = diag(3)) 
  }
}))

mat_1 <- scale(mat_1, center = T, scale = F)
mat_2 <- scale(mat_2, center = T, scale = F)
svd_1 <- svd(mat_1)
svd_2 <- svd(mat_2)

p_1 <- 10; p_2 <- 10
svd_v_1 <- tiltedCCA::generate_random_orthogonal(p_1, 3)
svd_v_2 <- tiltedCCA::generate_random_orthogonal(p_2, 3)

mat_1 <- tcrossprod(svd_1$u %*% diag(svd_1$d), svd_v_1)
mat_2 <- tcrossprod(svd_2$u %*% diag(svd_2$d), svd_v_2)

clustering_1 <- factor(stats::kmeans(mat_1, centers = 3)$cluster)
clustering_2 <- factor(stats::kmeans(mat_2, centers = 2)$cluster)

rownames(mat_1) <- paste0("n", 1:nrow(mat_1))
rownames(mat_2) <- paste0("n", 1:nrow(mat_2))
colnames(mat_1) <- paste0("g", 1:ncol(mat_1))
colnames(mat_2) <- paste0("p", 1:ncol(mat_2))

plot_idx <- sample(1:nrow(mat_1))
orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
blue_col <- rgb(129, 139, 191, maxColorValue = 255)
col_vec <- c(purple_col, orange_col, blue_col)

################################
# Step 2: Apply Tilted-CCA
################################

set.seed(10)
multiSVD_obj <- tiltedCCA::create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                           dims_1 = 1:3, dims_2 = 1:3,
                                           center_1 = F, center_2 = F,
                                           normalize_row = F,
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
                                        num_neigh = 80,
                                        bool_cosine = F,
                                        bool_intersect = F,
                                        min_deg = 0)
multiSVD_obj <- tiltedCCA::tiltedCCA(input_obj = multiSVD_obj,
                                     fix_tilt_perc = 0)
multiSVD_obj <- tiltedCCA::fine_tuning(input_obj = multiSVD_obj,
                                       verbose = 0)
multiSVD_obj <- tiltedCCA::tiltedCCA_decomposition(multiSVD_obj)

################################
# Step 3: Plot the data
################################

svd_func <- function(mat){
  svd_res <- svd(mat)
  svd_res$u[,1:2] %*% diag(svd_res$d[1:2])
}

png("simulation4_data-a.png", height = 600, width = 600, res = 300, units = "px")
par(mfrow = c(1,1), mar = c(2,2,0.5,0.5))
tmp <- svd_func(mat_1)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
graphics.off()

png("simulation4_data-b.png", height = 600, width = 600, res = 300, units = "px")
par(mfrow = c(1,1), mar = c(2,2,0.5,0.5))
tmp <- svd_func(mat_2)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
graphics.off()

################################
# Step 4: Plot Tilted-CCA's result
################################

names(multiSVD_obj)

height <- 1000
cex <- 1.3; cex.main <- 1.5
png("simulation4_tcca.png", height = height, width = 3000/1200*height, res = 300, units = "px")
par(mfrow = c(1,3))
tmp <- svd_func(multiSVD_obj$tcca_obj$common_score)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "Common embedding",
     xlab = "Common's dim. 1", ylab = "Common's dim. 2",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T,
     cex.lab = cex, cex = cex, cex.axis = cex, cex.main = cex.main)
tmp <- svd_func(multiSVD_obj$tcca_obj$distinct_score_1)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "Modality 1's distinct embedding",
     xlab = "Distinct-1's dim. 1", ylab = "Distinct-1's dim. 2",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T,
     cex.lab = cex, cex = cex, cex.axis = cex, cex.main = cex.main)
tmp <- svd_func(multiSVD_obj$tcca_obj$distinct_score_2)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "Modality 2's distinct embedding",
     xlab = "Distinct-2's dim. 1", ylab = "Distinct-2's dim. 2",
     pch = 16, col = true_cluster[plot_idx], asp = T,
     cex.lab = cex, cex = cex, cex.axis = cex, cex.main = cex.main)
graphics.off()

png("simulation4_tcca-a.png", height = 600, width = 600, res = 300, units = "px")
par(mfrow = c(1,1), mar = c(2,2,0.5,0.5))
tmp <- svd_func(multiSVD_obj$tcca_obj$common_score)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
graphics.off()


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

png("simulation4_consensuspca.png", height = 600, width = 600, res = 300, units = "px")
par(mfrow = c(1,1), mar = c(2,2,0.5,0.5))
tmp <- svd_func(consensus_pca$dimred_consensus)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
graphics.off()
