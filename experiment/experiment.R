rm(list=ls())

generate_custom_data <- function(seed = 10){
  set.seed(seed)
  n_clust <- 100
  B_mat <- matrix(c(0.9, 0.2, 0.1, 
                    0.2, 0.9, 0.1,
                    0.1, 0.1, 0.5), 3, 3, byrow = T)
  K <- ncol(B_mat); rho <- 1
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec); true_membership_vec <- membership_vec
  svd_u_1 <- generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)
  svd_u_2 <- generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)
  
  set.seed(seed)
  p_1 <- 20; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
  svd_v_1 <- generate_random_orthogonal(p_1, K-1)
  svd_v_2 <- generate_random_orthogonal(p_2, K-1)
  
  set.seed(seed)
  dat <- generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2)
  dcca_obj <- dcca_factor(dat$mat_1, dat$mat_2, dims_1 = 1:(K-1), dims_2 = 1:(K-1), 
                          verbose = F)
  membership_vec <- as.factor(membership_vec)
  set.seed(seed)
  list_g_1 <- construct_frnn(dcca_obj, nn = 5, membership_vec = membership_vec, 
                             data_1 = T, data_2 = F,
                             verbose = F, bool_matrix = T)
  set.seed(seed)
  list_g_2 <- construct_frnn(dcca_obj, nn = 5, membership_vec = membership_vec, 
                             data_1 = T, data_2 = F,
                             verbose = F, bool_matrix = T)
  
  set.seed(seed)
  c_g <- combine_frnn(dcca_obj, list_g_1$c_g, list_g_2$c_g, nn = 5)
  
  list(c_g = c_g, d_g_1 = list_g_1$d_g, d_g_2 = list_g_2$d_g)
}

#######

set.seed(10)
dat <- generate_custom_data()
c_g = dat$c_g
d_g_1 = dat$d_g_1
d_g_2 = dat$d_g_2
idx = 1:10
common_enrich = F
distinct_enrich_1 = T
distinct_enrich_2 = T
max_tries = 10

bool_vec <- c(common_enrich, distinct_enrich_1, distinct_enrich_2)
bool_idx <- which(bool_vec)
graph_list <- list(c_g, d_g_1, d_g_2)

neigh_list <- lapply(bool_idx, function(k){
  neigh_vec <- unlist(lapply(idx, function(i){.nonzero_col(graph_list[[k]], i, bool_value = F)}))
  neigh_vec <- setdiff(neigh_vec, idx)
})

all_idx <- sort(.common_intersection(neigh_list))
df <- sapply(1:length(neigh_list), function(k){
  table(intersect(neigh_list[[k]], all_idx))
})
min_neigh <- apply(df, 1, min)
