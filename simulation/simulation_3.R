rm(list=ls())
library(tiltedCCA)
set.seed(10)

################################
# Step 1: Generate the data, where Modality 1 is informative, but not Modality 2
################################

n_each <- 100
endpoint_list1 <- list(
  traj1 = matrix(c(0,0,0,
                   1,0,0,
                   1,0,3), 3, 3, byrow = T),
  traj2 = matrix(c(0,0,0,
                   3,0,0,
                   6,0,3), 3, 3, byrow = T),
  traj3 = matrix(c(0,0,0,
                   3,0,0,
                   6,0,-3), 3, 3, byrow = T)
)
endpoint_list2 <- list(
  traj1 = matrix(c(0,0,0,
                   1,0,0,
                   1,3,0), 3, 3, byrow = T),
  traj2 = matrix(c(0,0,0,
                   3,0,0,
                   6,0,0), 3, 3, byrow = T),
  traj3 = matrix(c(0,0,0,
                   3,0,0,
                   6,0,0), 3, 3, byrow = T)
)

length_func <- function(lis){
  lapply(lis, function(mat){
    p <- nrow(mat)
    vec <- sapply(2:p, function(j){
      sqrt(sum((mat[j,]-mat[j-1,])^2))
    })
    c(0, cumsum(vec/sum(vec)))
  })
}

segment_length_list1 <- length_func(endpoint_list1)
segment_length_list2 <- length_func(endpoint_list2)

sample_between_two_points <- function(point1, point2){
  rng <- stats::runif(1)
  rng*point1 + (1-rng)*point2
}

K <- ncol(endpoint_list1[[1]])
true_cluster <- rep(1:K, each = n_each)

generate_matrix <- function(endpoint_list,
                            K, 
                            n_each, 
                            segment_length_list,
                            sd){
  mat <- do.call(rbind, lapply(1:K, function(k){
    t(sapply(1:n_each, function(i){
      rng <- stats::runif(1)
      idx <- min(which(segment_length_list[[k]] > rng))
      sample_between_two_points(
        point1 = endpoint_list[[k]][idx-1,],
        point2 = endpoint_list[[k]][idx,]
      )
    }))
  }))
  mat <- mat + matrix(rnorm(prod(dim(mat)), mean = 0, sd = sd), 
                      ncol = ncol(mat), nrow = nrow(mat))
  mat
}

mat_1 <- generate_matrix(endpoint_list = endpoint_list1,
                         K = K, n_each = n_each,
                         segment_length_list = segment_length_list1,
                         sd = 0.1)
mat_2 <- generate_matrix(endpoint_list = endpoint_list2,
                         K = K, n_each = n_each,
                         segment_length_list = segment_length_list2,
                         sd = 0.1)

mat_1 <- scale(mat_1, center = T, scale = F)
mat_2 <- scale(mat_2, center = T, scale = F)
svd_1 <- svd(mat_1)
svd_2 <- svd(mat_2)

p_1 <- 10; p_2 <- 10
svd_v_1 <- tiltedCCA::generate_random_orthogonal(p_1, 3)
svd_v_2 <- tiltedCCA::generate_random_orthogonal(p_2, 3)

mat_1 <- tcrossprod(svd_1$u %*% diag(svd_1$d), svd_v_1)
mat_2 <- tcrossprod(svd_2$u %*% diag(svd_2$d), svd_v_2)

rownames(mat_1) <- paste0("n", 1:nrow(mat_1))
rownames(mat_2) <- paste0("n", 1:nrow(mat_2))
colnames(mat_1) <- paste0("g", 1:ncol(mat_1))
colnames(mat_2) <- paste0("p", 1:ncol(mat_2))

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
                                          large_clustering_1 = NULL, 
                                          large_clustering_2 = NULL,
                                          num_metacells = NULL)
multiSVD_obj <- tiltedCCA::compute_snns(input_obj = multiSVD_obj,
                                        latent_k = 2,
                                        num_neigh = 15,
                                        bool_cosine = F,
                                        bool_intersect = F,
                                        min_deg = 0)
multiSVD_obj <- tiltedCCA::tiltedCCA(input_obj = multiSVD_obj,
                                     enforce_boundary = F,
                                     fix_tilt_perc = 0)
multiSVD_obj <- tiltedCCA::fine_tuning(input_obj = multiSVD_obj,
                                       verbose = 0)
multiSVD_obj <- tiltedCCA::tiltedCCA_decomposition(multiSVD_obj)

################################
# Step 3: Plot the data
################################

# png("simulation3_data.png", height = 1200, width = 2000, res = 300, units = "px")
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
# image(t(multiSVD_obj$cca_obj$score_1))
# multiSVD_obj$tcca_obj$tilt_perc

svd_func <- function(mat){
  svd_res <- svd(mat)
  svd_res$u[,1:2] %*% diag(svd_res$d[1:2])
}

height <- 1000
cex <- 1.3; cex.main <- 1.5
# png("simulation3_tcca.png", height = height, width = 3000/1200*height, res = 300, units = "px")
par(mfrow = c(1,3))
tmp <- svd_func(cbind(multiSVD_obj$common_mat_1, multiSVD_obj$common_mat_2))
plot(tmp[,1], tmp[,2],
     main = "Common embedding",
     xlab = "Common's dim. 1", ylab = "Common's dim. 2",
     pch = 16, col = true_cluster, asp = T,
     cex.lab = cex, cex = cex, cex.axis = cex, cex.main = cex.main)
tmp <- svd_func(multiSVD_obj$distinct_mat_1)
plot(tmp[,1], tmp[,2],
     main = "Modality 1's distinct embedding",
     xlab = "Distinct-1's dim. 1", ylab = "Distinct-1's dim. 2",
     pch = 16, col = true_cluster, asp = T,
     cex.lab = cex, cex = cex, cex.axis = cex, cex.main = cex.main)
tmp <- svd_func(multiSVD_obj$distinct_mat_2)
plot(tmp[,1], tmp[,2],
     main = "Modality 2's distinct embedding",
     xlab = "Distinct-2's dim. 1", ylab = "Distinct-2's dim. 2",
     pch = 16, col = true_cluster, asp = T,
     cex.lab = cex, cex = cex, cex.axis = cex, cex.main = cex.main)
# graphics.off()

################################
# Step 5: For comparison, plot Consensus PCA
################################

set.seed(10)
consensus_pca <- tiltedCCA::consensus_pca(mat_1 = mat_1, mat_2 = mat_2,
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

png("simulation3_consensuspca.png", height = 1200, width = 1000, res = 300, units = "px")
par(mfrow = c(1,1))
plot(consensus_pca$dimred_consensus[,1], consensus_pca$dimred_consensus[,2],
     main = "Consensus PCA embedding",
     xlab = "Consensus PCA's dim. 1", ylab = "Consensus PCA's dim. 2",
     pch = 16, col = true_cluster, asp = T)
graphics.off()
