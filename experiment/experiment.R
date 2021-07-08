set.seed(10)
n_clust <- 100
B_mat <- matrix(c(0.9, 0.2, 0.1, 
                  0.2, 0.9, 0.1,
                  0.1, 0.1, 0.5), 3, 3, byrow = T)
K <- ncol(B_mat); rho <- 1
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
n <- length(membership_vec); true_membership_vec <- membership_vec
svd_u_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)
svd_u_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)

set.seed(10)
p_1 <- 20; p_2 <- 40
svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
svd_v_1 <- generate_random_orthogonal(p_1, K-1)
svd_v_2 <- generate_random_orthogonal(p_2, K-1)

set.seed(10)
dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2)
dcca_obj <- dcca_factor(dat$mat_1, dat$mat_2, dims_1 = 1:(K-1), dims_2 = 1:(K-1), 
                        verbose = F)
# res <- construct_frnn(dcca_obj, nn = 5, membership_vec = membership_vec, 
#                       verbose = F, bool_matrix = T)

###############3

membership_vec <- as.factor(membership_vec)
obj <- dcca_obj
nn = 5
data_1 = T
data_2 = F
max_subsample_frnn = nrow(obj$common_score)
frnn_approx = 0
radius_quantile = 0.9
bool_matrix = T
include_diag = T
verbose = T

embedding <- .prepare_embeddings(obj, data_1 = data_1, data_2 = data_2, 
                                 add_noise = F)
n <- nrow(embedding[[1]])
cell_subidx <- .construct_celltype_subsample(membership_vec, max_subsample_frnn)
if(length(cell_subidx) < n) {
  membership_vec <- membership_vec[cell_subidx]
}


