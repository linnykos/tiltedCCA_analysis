rm(list=ls())
n_each <- 100
true_membership_vec <- rep(1:3, each = n_each)
mat_1 <- do.call(rbind, lapply(1:3, function(i){
  if(i %in% c(1,2)){
    MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
  } else {
    MASS::mvrnorm(n = n_each, mu = c(9,0), Sigma = diag(2)) 
  }
}))

mat_2 <- do.call(rbind, lapply(1:3, function(i){
  if(i %in% c(1,3)){
    MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
  } else {
    MASS::mvrnorm(n = n_each, mu = c(9,0), Sigma = diag(2)) 
  }
}))

mat_1 <- scale(mat_1, center = T, scale = F)
mat_2 <- scale(mat_2, center = T, scale = F)
svd_1 <- svd(mat_1)
svd_2 <- svd(mat_2)

p_1 <- 40; p_2 <- 40
svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)

mat_1 <- tcrossprod(multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
mat_1 <- mat_1 + matrix(rnorm(prod(dim(mat_1)), sd = 1), ncol = ncol(mat_1), nrow = nrow(mat_1))
mat_2 <- tcrossprod(multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
mat_2 <- mat_2 + matrix(rnorm(prod(dim(mat_2)), sd = 1), ncol = ncol(mat_2), nrow = nrow(mat_2))

rownames(mat_1) <- paste0("c", 1:nrow(mat_1))
rownames(mat_2) <- paste0("c", 1:nrow(mat_2))
colnames(mat_1) <- paste0("g", 1:ncol(mat_1))
colnames(mat_2) <- paste0("p", 1:ncol(mat_2))

############

dims_1 = 1:2
dims_2 = 1:2
center_1 = F
center_2 = F
scale_1 = F
scale_2 = F
cell_max = nrow(mat_1)
clustering_resolution = 1
discretization_gridsize = 9
fix_tilt_perc = F
form_meta_matrix = F
metacell_clustering = as.factor(true_membership_vec)
num_neigh = min(30, round(nrow(mat_1)/20))
verbose = T

rank_1 <- max(dims_1); rank_2 <- max(dims_2)
stopifnot(nrow(mat_1) == nrow(mat_2), 
          rank_1 <= min(dim(mat_1)), rank_2 <= min(dim(mat_2)))

n <- nrow(mat_1)

if(verbose) print(paste0(Sys.time(),": D-CCA: Starting matrix shrinkage"))
svd_1 <- .svd_truncated(mat_1, K = rank_1, symmetric = F, rescale = F, 
                        mean_vec = center_1, sd_vec = scale_1, K_full_rank = F)
svd_2 <- .svd_truncated(mat_2, K = rank_2, symmetric = F, rescale = F, 
                        mean_vec = center_2, sd_vec = scale_2, K_full_rank = F)

svd_1 <- .check_svd(svd_1, dims = dims_1)
svd_2 <- .check_svd(svd_2, dims = dims_2)

metacell_clustering <- form_metacells(svd_1, svd_2, 
                                      clustering_resolution = clustering_resolution,
                                      dims_1 = NA, dims_2 = NA,
                                      center_1 = T, center_2 = T,
                                      scale_1 = T, scale_2 = T,
                                      verbose = verbose)
cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)

stopifnot(cell_max > 10)
full_rank <- length(cca_res$obj_vec)
n <- nrow(svd_1$u)

tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
score_1 <- tmp$score_1; score_2 <- tmp$score_2
stopifnot(ncol(score_1) == length(svd_1$d), ncol(score_2) == length(svd_2$d),
          nrow(score_1) == nrow(score_2))

obj_vec <- cca_res$obj_vec

# compute the common scores
n <- nrow(score_1)
if(cell_max < n){
  n_idx <- sample(1:n, size = cell_max)
} else {
  n_idx <- 1:n
}

# [[note to self: use these n_idx somehow]]
# tmp <- .common_decomposition(discretization_gridsize = discretization_gridsize,
#                              fix_tilt_perc = fix_tilt_perc,
#                              metacell_clustering = metacell_clustering,
#                              num_neigh = num_neigh,
#                              score_1 = score_1,
#                              score_2 = score_2,
#                              svd_1 = svd_1, 
#                              svd_2 = svd_2)

#####################

rank_c <- min(ncol(score_1), ncol(score_2))
stopifnot(all(sapply(1:rank_c, function(k){
  val <- score_1[,k] %*% score_2[,k]; val >= 0 
}))) # ensures score matrices contain pair of acute vectors

basis_list <- lapply(1:rank_c, function(k){
  .representation_2d(score_1[,k], score_2[,k])
})

circle_list <- lapply(1:rank_c, function(k){
  vec1 <- basis_list[[k]]$rep1
  vec2 <- basis_list[[k]]$rep2
  .construct_circle(vec1, vec2)
})

# if(verbose) print(paste0(Sys.time(),": D-CCA : (Inner) Computing distinct percentage"))
# if(is.logical(fix_tilt_perc) && !fix_tilt_perc){
#   tmp <- .search_tilt_perc(
#     basis_list = basis_list,
#     circle_list = circle_list,
#     discretization_gridsize = discretization_gridsize,
#     metacell_clustering = metacell_clustering,
#     num_neigh = num_neigh,
#     score_1 = score_1,
#     score_2 = score_2,
#     svd_1 = svd_1,
#     svd_2 = svd_2
#   )
#   tilt_perc <- tmp$percentage
#   df_percentage <- tmp$df
# } 

r <- length(basis_list)

# handle corner case
tol <- 1e-6
if(sum(sapply(1:r, function(k){
  sum(abs(basis_list[[k]]$rep1 - basis_list[[k]]$rep2))
})) <= tol) return(0.5)

# initialize values
percentage_grid <- seq(0, 1, length.out = discretization_gridsize)
value_vec <- rep(NA, discretization_gridsize)
percentage <- percentage_grid[5]

r <- length(basis_list)

radian_vec <- sapply(1:r, function(k){
  .compute_radian(percentage, 
                  vec1 = basis_list[[k]]$rep1,
                  vec2 = basis_list[[k]]$rep2,
                  circle = circle_list[[k]])
})

common_representation <- sapply(1:r, function(k){
  .position_from_circle(circle_list[[k]], radian_vec[k])
})

common_score <- sapply(1:r, function(k){
  basis_list[[k]]$basis_mat %*% common_representation[,k]
})

common_mat <- .convert_common_score_to_mat(common_score,
                                           score_1,
                                           score_2,
                                           svd_1, 
                                           svd_2)

mat <- common_mat 
n <- nrow(mat)
snn_mat <- .form_snn_mat(bool_intersect = T,
                         mat, num_neigh = num_neigh)
image(as.matrix(snn_mat), asp = T)
density_mat <- .compute_density_matrix(as.matrix(snn_mat),
                                       metacell_clustering = metacell_clustering)
K <- length(levels(metacell_clustering))

tol <- 1e-3
quality_vec <- sapply(1:K, function(k){
  offdiag_mean <- mean(density_mat[k,-k])
  offdiag_sd <- sd(density_mat[k,-k])
  # print(paste0(k, ": ", round(density_mat[k,k],2)))
  # print(paste0(k, ": ", round(offdiag_mean+offdiag_sd,2)))
  
  density_mat[k,k]/max(c(offdiag_mean-offdiag_sd, tol))
})

mean(quality_vec) * sum(snn_mat)/n


