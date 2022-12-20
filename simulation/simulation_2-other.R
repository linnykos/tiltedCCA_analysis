rm(list=ls())
library(tiltedCCA)
set.seed(10)
source("../tiltedCCA_analysis/main/jive.R")

################################
# Step 1: Generate the data, where Modality 1 is informative, but not Modality 2
################################

n_each <- 100
true_cluster <- rep(1:5, each = n_each)
mat_1 <- do.call(rbind, lapply(1:5, function(i){
  if(i %in% c(1,2)){
    MASS::mvrnorm(n = n_each, mu = c(10,10,0), Sigma = diag(3))  
  } else if(i %in% c(3,4)) {
    MASS::mvrnorm(n = n_each, mu = c(0,0,0), Sigma = diag(3)) 
  } else {
    MASS::mvrnorm(n = n_each, mu = c(0,0,10), Sigma = diag(3)) 
  }
}))

mat_2 <- do.call(rbind, lapply(1:5, function(i){
  if(i %in% c(1,3)){
    MASS::mvrnorm(n = n_each, mu = c(10,10,0), Sigma = diag(3)) 
  } else if(i %in% c(2,4)) {
    MASS::mvrnorm(n = n_each,mu = c(0,0,0), Sigma = diag(3)) 
  } else {
    MASS::mvrnorm(n = n_each, mu = c(0,0,10), Sigma = diag(3)) 
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
clustering_2 <- factor(stats::kmeans(mat_2, centers = 3)$cluster)

rownames(mat_1) <- paste0("n", 1:nrow(mat_1))
rownames(mat_2) <- paste0("n", 1:nrow(mat_2))
colnames(mat_1) <- paste0("g", 1:ncol(mat_1))
colnames(mat_2) <- paste0("p", 1:ncol(mat_2))

plot_idx <- sample(1:nrow(mat_1))
orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
blue_col <- rgb(129, 139, 191, maxColorValue = 255)
black_col <- "black"
green_col <- rgb(70, 177, 70, maxColorValue = 255)
col_vec <- c(orange_col, purple_col, blue_col, black_col, green_col)

################################

svd_func <- function(mat){
  svd_res <- svd(mat)
  svd_res$u[,1:2] %*% diag(svd_res$d[1:2])
}

par(mfrow = c(1,2), mar = c(2,2,0.5,0.5))
tmp <- svd_func(mat_1)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
tmp <- svd_func(mat_2)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)


################################
# JIVE
################################

set.seed(10)
jive_results <- jive(mat_1 = mat_1, 
                     mat_2 = mat_2, 
                     common_r = 3,
                     r_1 = 3,
                     r_2 = 3,
                     return_full_prediction = T)

names(jive_results)

# png("simulation4_consensuspca.png", height = 600, width = 600, res = 300, units = "px")
par(mfrow = c(1,1), mar = c(2,2,0.5,0.5))
tmp <- svd_func(jive_results$embedding)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
# graphics.off()

# plot the reconstruction
par(mfrow = c(1,2), mar = c(2,2,0.5,0.5))
tmp <- svd_func(jive_results$pred_mat_1 + jive_results$a_mat_1)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
tmp <- svd_func(jive_results$pred_mat_2 + jive_results$a_mat_2)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)

# plot the individual embeddings
par(mfrow = c(1,2), mar = c(2,2,0.5,0.5))
tmp <- svd_func(jive_results$a_embedding_1)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
tmp <- svd_func(jive_results$a_embedding_2)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)

# plot the embedding of the individual matrices
par(mfrow = c(1,2), mar = c(2,2,0.5,0.5))
tmp <- svd_func(jive_results$a_mat_1)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
tmp <- svd_func(jive_results$a_mat_2)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)

# plot the embedding of the joint matrices
par(mfrow = c(1,2), mar = c(2,2,0.5,0.5))
tmp <- svd_func(jive_results$pred_mat_1)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
tmp <- svd_func(jive_results$pred_mat_2)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)

