rm(list=ls())
n <- 100; K <- 2
common_space <- scale(MASS::mvrnorm(n = n, mu = rep(0,K), Sigma = diag(K)), center = T, scale = F)

p1 <- 5; p2 <- 10
transform_mat_1 <- matrix(stats::runif(K*p1, min = -1, max = 1), nrow = K, ncol = p1)
transform_mat_2 <- matrix(stats::runif(K*p2, min = -1, max = 1), nrow = K, ncol = p2)

mat_1 <- common_space %*% transform_mat_1; mat_2 <- common_space %*% transform_mat_2 

# dcca_res <- dcca_factor(mat_1, mat_2, dims_1 = 1:K, dims_2 = 1:K, 
#                         verbose = F, center_1 = F, center_2 = F, 
#                         scale_1 = F, scale_2 = F)
###########################
dims_1 = 1:K
dims_2 = 1:K
verbose = F
center_1 = F
center_2 = F
scale_1 = F
scale_2 = F
cell_max = nrow(mat_1)
discretization_gridsize = 9
enforce_boundary = is.factor(metacell_clustering_1)
fix_tilt_perc = F
form_meta_matrix = F
metacell_clustering_1 = NA
metacell_clustering_2 = NA
num_neigh = min(30, round(nrow(mat_1)/20))
verbose = T

stopifnot((all(is.na(metacell_clustering_1)) & all(is.na(metacell_clustering_2))) ||
            (is.list(metacell_clustering_1) & is.list(metacell_clustering_2)) ||
            (is.factor(metacell_clustering_1) & is.factor(metacell_clustering_2)))
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

if(all(is.na(metacell_clustering_1)) & all(is.na(metacell_clustering_2))){
  tmp <- .form_snns(num_neigh = num_neigh, svd_1 = svd_1, svd_2 = svd_2)
  metacell_clustering_1 <- tmp$metacell_clustering_1
  metacell_clustering_2 <- tmp$metacell_clustering_2
}

cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)

res <- .dcca_common_score(cca_res = cca_res, 
                          cell_max = cell_max,
                          check_alignment = all(!is.na(metacell_clustering_1)), 
                          discretization_gridsize = discretization_gridsize,
                          enforce_boundary = T,
                          fix_tilt_perc = T, 
                          metacell_clustering_1 = metacell_clustering_1,
                          metacell_clustering_2 = metacell_clustering_2,
                          num_neigh = num_neigh,
                          svd_1 = svd_1, 
                          svd_2 = svd_2, 
                          verbose = verbose, msg = "")

############################3

stopifnot(cell_max > 10)
full_rank <- length(cca_res$obj_vec)
n <- nrow(svd_1$u)

if(verbose) print(paste0(Sys.time(),": D-CCA", msg, ": Computing unnormalized scores"))
tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
score_1 <- tmp$score_1; score_2 <- tmp$score_2
stopifnot(ncol(score_1) == length(svd_1$d), ncol(score_2) == length(svd_2$d),
          nrow(score_1) == nrow(score_2))

n <- nrow(score_1)
if(cell_max < n){
  n_idx <- sample(1:n, size = cell_max)
} else {
  n_idx <- 1:n
}
tmp <- .common_decomposition(discretization_gridsize = discretization_gridsize,
                             enforce_boundary = enforce_boundary,
                             fix_tilt_perc = fix_tilt_perc,
                             metacell_clustering_1 = metacell_clustering_1,
                             metacell_clustering_2 = metacell_clustering_2,
                             n_idx = n_idx,
                             num_neigh = num_neigh,
                             score_1 = score_1,
                             score_2 = score_2,
                             svd_1 = svd_1, 
                             svd_2 = svd_2,
                             verbose = verbose)

#########################3

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



