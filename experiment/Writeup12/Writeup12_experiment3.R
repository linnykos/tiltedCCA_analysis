rm(list=ls())
library(Seurat)
source("../multiomicCCA_analysis/experiment/Writeup12/Writeup12_simulation_functions.R")

set.seed(10)
n_clust <- 100
high <- 0.9; low <- 0.05
B_mat1 <- matrix(c(high, high, low, 
                   high, high, low,
                   low, low, high), 3, 3)
B_mat2 <- matrix(c(high, low, low, 
                   low, high, high,
                   low, high, high), 3, 3)
K <- ncol(B_mat1)
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
n <- length(membership_vec); true_membership_vec <- membership_vec
score_1 <- generate_sbm_orthogonal(B_mat1, membership_vec)
score_2 <- generate_sbm_orthogonal(B_mat2, membership_vec)
#score_2 <- score_1[c((2*n_clust+1):(3*n_clust),1:(2*n_clust)),]

set.seed(10)
p_1 <- 20; p_2 <- 20
coef_mat1 <- matrix(stats::rnorm(K*p_1), K, p_1)
coef_mat2 <- matrix(stats::rnorm(K*p_1), K, p_1)

set.seed(10)
dat <- generate_data(score_1, score_2, coef_mat1, coef_mat2, 
                     noise_func = function(mat){matrix(stats::rnorm(prod(dim(mat)), mean = mat, sd = 1.5), nrow(mat), ncol(mat))})

set.seed(10)
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, apply_shrinkage = F, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

prep_list <- .prepare_umap_embedding(dcca_decomp)
zz1 <- .extract_matrix_helper(prep_list$common_score, prep_list$distinct_score_1,
                             prep_list$svd_list$e1, common_bool = T, distinct_bool = F)
zz2 <- .extract_matrix_helper(prep_list$common_score, prep_list$distinct_score_1,
                              prep_list$svd_list$e1, common_bool = F, distinct_bool = T)
zz <- rbind(zz1, zz2)
set.seed(10)
tmp <- Seurat::RunUMAP(zz, metric = "euclidean", verbose = F)@cell.embeddings
par(mfrow = c(1,2))
plot(tmp[1:n,1], tmp[1:n,2], col = true_membership_vec, pch = 16, main = "Common")
plot(tmp[(n+1):(2*n),1], tmp[(n+1):(2*n),2], col = true_membership_vec, pch = 16, main = "Distinct")

##########3


par(mfrow = c(1,3))
image(.rotate(zz))
plot(zz[,1])
image(.rotate(prep_list$common_score))

par(mfrow = c(1,2))
set.seed(10)
tmp <- Seurat::RunUMAP(zz, metric = "euclidean", verbose = F, min.dist = 0.3)@cell.embeddings
plot(tmp[,1], tmp[,2], col = true_membership_vec, pch = 16)
set.seed(10)
tmp <- Seurat::RunUMAP(zz, metric = "euclidean", verbose = F, a = 1, b = 1e-5)@cell.embeddings
plot(tmp[,1], tmp[,2], col = true_membership_vec, pch = 16)

#####

zz <- .extract_matrix_helper(prep_list$common_score, prep_list$distinct_score_1,
                             prep_list$svd_list$e1, common_bool = F, distinct_bool = T)
par(mfrow = c(1,3))
image(.rotate(zz))
plot(zz[,1])
image(.rotate(prep_list$distinct_score_1))

########################

common_score = prep_list$common_score
distinct_score = prep_list$distinct_score_1
svd_e = prep_list$svd_list$e1
common_bool = T
distinct_bool = F

full_mat <- .mult_mat_vec(svd_e$u, svd_e$d)
center_vec <- apply(full_mat, 2, mean)
full_mat <- sapply(1:ncol(full_mat), function(k){
  full_mat[,k] - center_vec[k]
})
l2_vec <- apply(full_mat, 1, function(x){.l2norm(x)})
full_mat2 <- .mult_vec_mat(1/l2_vec, full_mat)

if(ncol(common_score) < ncol(distinct_score)){
  common_score <- rbind(common_score, matrix(0, nrow = nrow(common_score), ncol = ncol(distinct_score) - ncol(common_score)))
}
canonical_score <- common_score+distinct_score

if(common_bool){ 
  tmp <- tcrossprod(common_score, canonical_score) %*% full_mat
} else { 
  tmp <- tcrossprod(distinct_score, canonical_score) %*% full_mat
}

# center variables
tmp <- sapply(1:ncol(tmp), function(k){
  tmp[,k] - center_vec[k]
})

tmp <- .mult_vec_mat(1/l2_vec, tmp)

