
# high distinct 1, low distinct 2 with 3 clusters
set.seed(10)
n_clust <- 100
B_mat1 <- matrix(c(0.9, 0, 0, 
                   0, 0.9, 0,
                   0, 0, 0.9), 3, 3, byrow = T)
B_mat2 <- matrix(c(0.9, 0.85, 0, 
                   0.85, 0.9, 0,
                   0, 0, 1), 3, 3, byrow = T)
K <- ncol(B_mat1)
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
n <- length(membership_vec); true_membership_vec <- membership_vec
svd_u_1 <- generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)
svd_u_2 <- generate_sbm_orthogonal(B_mat2, membership_vec, centered = T)

set.seed(10)
p_1 <- 20; p_2 <- 40
svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
svd_v_1 <- generate_random_orthogonal(p_1, K-1)
svd_v_2 <- generate_random_orthogonal(p_2, K-1)

# set.seed(10)
# dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2)

############

num_neigh = max(round(nrow(svd_u_1)/20), 40)
noise_val = 1

n <- nrow(svd_u_1); p1 <- nrow(svd_v_1); p2 <- nrow(svd_v_2)
agg_mat <- crossprod(svd_u_1, svd_u_2)
cca_svd <- svd(agg_mat)
score_1 <- svd_u_1 %*% cca_svd$u
score_2 <- svd_u_2 %*% cca_svd$v

mat_1 <- tcrossprod(.mult_mat_vec(svd_u_1, svd_d_1), svd_v_1)
mat_2 <- tcrossprod(.mult_mat_vec(svd_u_2, svd_d_2), svd_v_2)
nn_1 <- RANN::nn2(.mult_mat_vec(svd_u_1, svd_d_1), k = num_neigh)$nn.idx
nn_2 <- RANN::nn2(.mult_mat_vec(svd_u_2, svd_d_2), k = num_neigh)$nn.idx

# res <- .common_decomposition(score_1, score_2, nn_1, nn_2)
#############

fix_common_perc = F
tol = 1e-6

stopifnot((!any(is.na(nn_1)) & !any(is.na(nn_2)) & !fix_common_perc) | 
            (all(is.na(nn_1)) & all(is.na(nn_2)) & fix_common_perc))

rank_c <- min(ncol(score_1), ncol(score_2))
stopifnot(all(sapply(1:rank_c, function(k){
  val <- score_1[,k] %*% score_2[,k]; val >= 0 
}))) # ensures score matrices contain pair of acute vectors

basis_list <- lapply(1:rank_c, function(k){
  .representation_2d(score_1[,k], score_2[,k])
})

.latent_common_perc(score_1[,1], score_2[,1], nn_1, nn_2)
plot(score_1[,1], score_2[,1], asp = T)

# if(fix_common_perc){
#   common_perc <- rep(0.5, rank_c)
# } else {
#   common_perc <- sapply(1:rank_c, function(k){
#     if(sum(abs(basis_list[[k]]$rep1 - basis_list[[k]]$rep2)) < tol){
#       .5
#     } else {
#       .latent_common_perc(score_1[,k], score_2[,k], nn_1, nn_2)
#     }
#   })
# }

#33333333#3

score_vec_1 <- score_1[,1]; score_vec_2 <- score_2[,1]
stopifnot(length(score_vec_1) == length(score_vec_2))

n <- length(score_vec_1)
mode1_common_perc <- sapply(1:n, function(i){
  print(i)
  val_1in1 <- stats::sd(score_vec_1[i] - score_vec_1[nn_1[i,]])
  val_1in2 <- stats::sd(score_vec_1[i] - score_vec_1[nn_2[i,]])
  val_2in1 <- stats::sd(score_vec_2[i] - score_vec_2[nn_1[i,]])
  val_2in2 <- stats::sd(score_vec_2[i] - score_vec_2[nn_2[i,]])
  
  .sigmoid_ratio(max(val_1in2/val_1in1, 1), max(val_2in1/val_2in2, 1))
})
