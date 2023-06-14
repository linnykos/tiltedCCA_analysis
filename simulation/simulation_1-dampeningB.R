rm(list=ls())
library(tiltedCCA)
library(mclust)
set.seed(10)

################################
# Step 1: Generate the data, where Modality 1 is informative, but not Modality 2
################################

simulation_function <- function(separation = 1){
  n_each <- 100
  B_mat1 <- matrix(c(0.9, 0.1, 0.1,
                     0.1, 0.9, 0.1,
                     0.1, 0.1, 0.9), 3, 3, byrow = T)
  K <- ncol(B_mat1)
  membership_vec <- rep(1:3, each = n_each)
  true_cluster <- membership_vec
  n <- length(membership_vec); true_membership_vec <- membership_vec
  svd_u_1 <- tiltedCCA::generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)[,1:2]
  svd_u_2a <- tiltedCCA::generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)[,1:2]
  svd_u_2b <- tiltedCCA::generate_random_orthogonal(n, 2, centered = T)
  svd_u_2 <- separation*svd_u_2a + (1-separation)*svd_u_2b
  svd_u_2[,1] <- svd_u_2[,1]/sqrt(sum(svd_u_2[,1]^2))
  svd_u_2[,2] <- svd_u_2[,2]/sqrt(sum(svd_u_2[,2]^2))
    
  p_1 <- 10; p_2 <- 10
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
  svd_v_1 <- tiltedCCA::generate_random_orthogonal(p_1, 2)
  svd_v_2 <- tiltedCCA::generate_random_orthogonal(p_2, 2)
  
  mat_1 <- tcrossprod(svd_u_1 %*% diag(svd_d_1), svd_v_1)
  mat_2 <- tcrossprod(svd_u_2 %*% diag(svd_d_2), svd_v_2)
  
  mat_1 <- scale(mat_1, center = T, scale = T)
  mat_2 <- scale(mat_2, center = T, scale = T)
  svd_2 <- svd(mat_2)
  
  # number of clusters for mat_2
  mclust_res <- mclust::Mclust(svd_2$u %*% diag(svd_2$d), 
                               G = 1:3,
                               modelNames = "VVV",
                               verbose = F)
  
  clustering_1 <- factor(stats::kmeans(mat_1, centers = 3)$cluster)
  clustering_2 <- factor(stats::kmeans(mat_2, centers = mclust_res$G)$cluster)
  
  rownames(mat_1) <- paste0("n", 1:nrow(mat_1))
  rownames(mat_2) <- paste0("n", 1:nrow(mat_2))
  colnames(mat_1) <- paste0("g", 1:ncol(mat_1))
  colnames(mat_2) <- paste0("p", 1:ncol(mat_2))
  
  list(clustering_1 = clustering_1,
       clustering_2 = clustering_2,
       mat_1 = mat_1,
       mat_2 = mat_2,
       true_cluster = true_cluster)
}

svd_func <- function(mat){
  svd_res <- svd(mat)
  svd_res$u[,1:2] %*% diag(svd_res$d[1:2])
}

plot_idx <- sample(1:300) # 300 is equal to nrow(mat_1)
orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
blue_col <- rgb(129, 139, 191, maxColorValue = 255)
col_vec <- c(purple_col, orange_col, blue_col)

################################
# Step 2: Apply Tilted-CCA for sequence of separations: 100 trials for 5 different values
################################

trials <- 50

separation_vec <- seq(sqrt(0), sqrt(1), length.out=10)
len <- length(separation_vec)
result_list <- lapply(1:len, function(i){
  list(clustering_1 = NA,
       clustering_2 = NA,
       mat_1 = NA,
       mat_2 = NA,
       example_res = NA, 
       true_cluster = NA,
       clusters = rep(NA, trials))
})
names(result_list) <- paste0("separation_", separation_vec)

for(separation in separation_vec){
  print("=====")
  print(paste0("Working on separation ", round(separation,2)))
  for(trial in 1:trials){
    set.seed(10*trial)
    if(trial %% floor(trials/10) == 0) cat('*')
    
    # simulate data
    res <- simulation_function(separation = separation)
    
    set.seed(10)
    multiSVD_obj <- tiltedCCA::create_multiSVD(mat_1 = res$mat_1, mat_2 = res$mat_2,
                                               dims_1 = 1:2, dims_2 = 1:2,
                                               center_1 = T, center_2 = T,
                                               normalize_row = F,
                                               normalize_singular_value = F,
                                               recenter_1 = F, recenter_2 = F,
                                               rescale_1 = F, rescale_2 = F,
                                               scale_1 = T, scale_2 = T)
    multiSVD_obj <- tiltedCCA::form_metacells(input_obj = multiSVD_obj,
                                              large_clustering_1 = res$clustering_1, 
                                              large_clustering_2 = res$clustering_2,
                                              num_metacells = NULL)
    multiSVD_obj <- tiltedCCA::compute_snns(input_obj = multiSVD_obj,
                                            latent_k = 2,
                                            num_neigh = 80,
                                            bool_cosine = F,
                                            bool_intersect = F,
                                            min_deg = 0)
    multiSVD_obj <- tiltedCCA::tiltedCCA(input_obj = multiSVD_obj,
                                         fix_tilt_perc = F)
    multiSVD_obj <- tiltedCCA::fine_tuning(input_obj = multiSVD_obj,
                                           verbose = 0)
    multiSVD_obj <- tiltedCCA::tiltedCCA_decomposition(multiSVD_obj)
    
    if(trial == 1){
      result_list[[paste0("separation_", separation)]]$clustering_1 <- res$clustering_1
      result_list[[paste0("separation_", separation)]]$clustering_2 <- res$clustering_2
      result_list[[paste0("separation_", separation)]]$mat_1 <- res$mat_1
      result_list[[paste0("separation_", separation)]]$mat_2 <- res$mat_2
      result_list[[paste0("separation_", separation)]]$true_cluster <- res$true_cluster
      result_list[[paste0("separation_", separation)]]$example_res <- multiSVD_obj
    }
    
    tmp <- svd_func(multiSVD_obj$tcca_obj$common_score)
    
    mclust_res <- mclust::Mclust(tmp, 
                                 G = 1:3,
                                 modelNames = "VVV",
                                 verbose = F)
    result_list[[paste0("separation_", separation)]]$clusters[trial] <- mclust_res$G
  }
  
  print("\n")
  print(table(result_list[[paste0("separation_", separation)]]$clusters))
}

###########

for(i in 1:length(result_list)){
  print(table(result_list[[i]]$clusters))
}

###########

mat_1 <- result_list[[1]]$mat_1
mat_2 <- result_list[[1]]$mat_2
true_cluster <- result_list[[1]]$true_cluster
par(mfrow = c(1,2))
tmp <- svd_func(mat_1)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "Data 1",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
tmp <- svd_func(mat_2)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "Data 2",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)

###########

cex <- 1.3; cex.main <- 1.5
multiSVD_obj <- result_list[[6]]$example_res
true_cluster <- result_list[[6]]$true_cluster
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
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T,
     cex.lab = cex, cex = cex, cex.axis = cex, cex.main = cex.main)




