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

########################################

input_obj = multiSVD_obj
discretization_gridsize = 21
enforce_boundary = F
fix_tilt_perc = 1
verbose = 0

if(verbose >= 1) print(paste0(Sys.time(),": Tilted-CCA: Gathering relevant objects"))
input_obj <- .set_defaultAssay(input_obj, assay = 1)
svd_1 <- .get_SVD(input_obj)
input_obj <- .set_defaultAssay(input_obj, assay = 2)
svd_2 <- .get_SVD(input_obj)

n <- nrow(svd_1$u)
metacell_clustering_list <- .get_metacell(input_obj,
                                          resolution = "cell", 
                                          type = "list", 
                                          what = "metacell_clustering")
if(!all(is.null(metacell_clustering_list))){
  averaging_mat <- .generate_averaging_matrix(metacell_clustering_list = metacell_clustering_list,
                                              n = n)
} else {
  averaging_mat <- NULL
}

target_dimred <- .get_Laplacian(input_obj, bool_common = T)
param <- .get_param(input_obj)
snn_bool_cosine <- param$snn_bool_cosine
snn_bool_intersect <- param$snn_bool_intersect
snn_k <- param$snn_latent_k
snn_min_deg <- param$snn_min_deg
snn_num_neigh <- param$snn_num_neigh

if(verbose >= 1) print(paste0(Sys.time(),": Tilted-CCA: Computing CCA"))
cca_res <- .cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)

# res <- .tiltedCCA_common_score(averaging_mat = averaging_mat,
#                                cca_res = cca_res, 
#                                discretization_gridsize = discretization_gridsize,
#                                enforce_boundary = enforce_boundary,
#                                fix_tilt_perc = fix_tilt_perc, 
#                                snn_bool_cosine = snn_bool_cosine,
#                                snn_bool_intersect = snn_bool_intersect,
#                                snn_k = snn_k,
#                                snn_min_deg = snn_min_deg,
#                                snn_num_neigh = snn_num_neigh,
#                                svd_1 = svd_1, 
#                                svd_2 = svd_2, 
#                                target_dimred = target_dimred,
#                                verbose = verbose)
######################

full_rank <- length(cca_res$obj_vec)
n <- nrow(svd_1$u)

if(verbose >= 1) print(paste0(Sys.time(),": Tilted-CCA: Computing unnormalized scores"))
tmp <- .compute_unnormalized_scores(svd_1, svd_2, cca_res)
score_1 <- tmp$score_1; score_2 <- tmp$score_2
stopifnot(ncol(score_1) == length(svd_1$d), ncol(score_2) == length(svd_2$d),
          nrow(score_1) == nrow(score_2))

if(verbose >= 1) print(paste0(Sys.time(),": Tilted-CCA: Computing common factors"))
obj_vec <- cca_res$obj_vec

# tmp <- .common_decomposition(averaging_mat = averaging_mat,
#                              discretization_gridsize = discretization_gridsize,
#                              enforce_boundary = enforce_boundary,
#                              fix_tilt_perc = fix_tilt_perc,
#                              score_1 = score_1,
#                              score_2 = score_2,
#                              snn_bool_cosine = snn_bool_cosine,
#                              snn_bool_intersect = snn_bool_intersect,
#                              snn_k = snn_k,
#                              snn_min_deg = snn_min_deg,
#                              snn_num_neigh = snn_num_neigh,
#                              svd_1 = svd_1, 
#                              svd_2 = svd_2,
#                              target_dimred = target_dimred,
#                              verbose = verbose)

######################

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

tilt_perc <- fix_tilt_perc; df_percentage <- NA

# tmp <- .evaluate_radian(
#   averaging_mat = averaging_mat,
#   basis_list = basis_list, 
#   circle_list = circle_list,
#   enforce_boundary = enforce_boundary,
#   percentage = tilt_perc,
#   return_common_score_basis = T,
#   score_1 = score_1,
#   score_2 = score_2,
#   snn_bool_cosine = snn_bool_cosine,
#   snn_bool_intersect = snn_bool_intersect,
#   snn_k = snn_k,
#   snn_min_deg = snn_min_deg,
#   snn_num_neigh = snn_num_neigh,
#   svd_1 = svd_1, 
#   svd_2 = svd_2,
#   target_dimred = target_dimred,
#   verbose = verbose
# )

#################

percentage = tilt_perc
return_common_score_basis = T

r <- length(basis_list)

radian_vec <- sapply(1:r, function(k){
  .compute_radian(circle = circle_list[[k]],
                  enforce_boundary = enforce_boundary,
                  percentage_val = 1, 
                  vec1 = basis_list[[k]]$rep1,
                  vec2 = basis_list[[k]]$rep2)
})

common_representation <- sapply(1:r, function(k){
  .position_from_circle(circle_list[[k]], radian_vec[k])
})

common_score <- sapply(1:r, function(k){
  basis_list[[k]]$basis_mat %*% common_representation[,k]
})
