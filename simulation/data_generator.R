simulation_all <- function(idx){
  if(idx == 1) { simulation1() 
  } else if(idx == 2){ simulation2() 
  } else if(idx == 3){ simulation3() 
  } else if(idx == 4){ simulation4() 
  } else if(idx == 5){ simulation5() 
  } else if(idx == 6){ simulation6() 
  } else {
    stop()
  }
}

# easy setting: low distinct 1, low distinct 2 with 3 clusters
simulation1 <- function(){
  n_clust <- 100
  B_mat <- matrix(c(0.9, 0.2, 0.1, 
                    0.2, 0.9, 0.1,
                    0.1, 0.1, 0.5), 3, 3, byrow = T)
  K <- ncol(B_mat); rho <- 1
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec); true_membership_vec <- membership_vec
  svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)
  svd_u_2 <- multiomicCCA::generate_sbm_orthogonal(rho*B_mat, membership_vec, centered = T)
  
  p_1 <- 20; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
  svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, K-1)
  svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, K-1)
  
  dat <- multiomicCCA::generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2)
  list(dat = dat, true_membership_vec = true_membership_vec)
}

# high distinct 1, low distinct 2 with 3 clusters
simulation2 <- function(){
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
  svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)
  svd_u_2 <- multiomicCCA::generate_sbm_orthogonal(B_mat2, membership_vec, centered = T)
  
  p_1 <- 20; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
  svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, K-1)
  svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, K-1)
  
  dat <- multiomicCCA::generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2, 
                                     noise_val = 0.1)
  list(dat = dat, true_membership_vec = true_membership_vec)
}

# even more dramatic high distinct 1, low distinct 2 with 3 clusters
simulation3 <- function(){
  n_clust <- 100
  high <- 0.9; low <- 0.05
  B_mat1 <- matrix(c(0.9, 0.1, 0.1,
                     0.1, 0.9, 0.1,
                     0.1, 0.1, 0.9), 3, 3, byrow = T)
  K <- ncol(B_mat1)
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec); true_membership_vec <- membership_vec
  svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)[,1:2]
  svd_u_2 <- multiomicCCA::generate_random_orthogonal(n, 2, centered = T)
  
  p_1 <- 20; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
  svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, 2)
  svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, 2)
  
  dat <- multiomicCCA::generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_1, 
                                     noise_val = 2)
  list(dat = dat, true_membership_vec = true_membership_vec)
}

# high distinct 1, high distinct 2 with 4 clusters. 
# distinct 1 separates (1) from (2) and 
# distinct 2 seperates (3) from (4) analogously
simulation4 <- function(){
  B_mat <- matrix(c(0.9, 0.4, 0.1,
                    0.4, 0.9, 0.1,
                    0.1, 0.1, 0.3), 3, 3, byrow = T)
  K <- ncol(B_mat); n_clust <- 100
  
  true_membership_vec <- rep(1:4, each = n_clust)
  n <- length(true_membership_vec)
  
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, 2*n_clust))
  svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(B_mat, membership_vec, centered = T)
  membership_vec <- c(rep(3, 2*n_clust), rep(1, n_clust), rep(2, n_clust))
  svd_u_2 <- multiomicCCA::generate_sbm_orthogonal(B_mat, membership_vec, centered = T)

  p_1 <- 20; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
  svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, K-1)
  svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, K-1)
  
  dat <- multiomicCCA::generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2)
  list(dat = dat, true_membership_vec = true_membership_vec)
}

# high distinct 1, high distinct 2 with 5 clusters
# the common loading separates (1234) from (5)
# in dataset 1, we can split black&red (12) from green&blue(34)
# in dataset 2, we can split black&green (13) from red&blue(24)
simulation5 <- function(){
  B_mat <- matrix(c(0.9, 0.4, 0.1, 
                    0.4, 0.9, 0.1, 
                    0.1, 0.1, 0.9), 3, 3, byrow = T)
  K <- ncol(B_mat); n_clust <- 100
  true_membership_vec <- rep(1:5, each = n_clust)
  n <- length(true_membership_vec)
  
  membership_vec <-  c(rep(1, 2*n_clust), rep(2, 2*n_clust), rep(3, n_clust))
  svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(B_mat, membership_vec, centered = T)
  membership_vec <-  c(rep(1, n_clust), rep(2, n_clust), rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  svd_u_2 <- multiomicCCA::generate_sbm_orthogonal(B_mat, membership_vec, centered = T)
  
  p_1 <- 40; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
  svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, K-1)
  svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, K-1)
  
  dat <- multiomicCCA::generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2,
                       noise_val = 2)
  list(dat = dat, true_membership_vec = true_membership_vec)
}

# same setting as above, but now w/ different number variables
simulation6 <- function(){
  B_mat <- matrix(c(0.9, 0.4, 0.1, 
                    0.4, 0.9, 0.1, 
                    0.1, 0.1, 0.9), 3, 3, byrow = T)
  K <- ncol(B_mat); n_clust <- 100
  true_membership_vec <- rep(1:5, each = n_clust)
  n <- length(true_membership_vec)
  
  membership_vec <-  c(rep(1, 2*n_clust), rep(2, 2*n_clust), rep(3, n_clust))
  svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(B_mat, membership_vec, centered = T)
  membership_vec <-  c(rep(1, n_clust), rep(2, n_clust), rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  svd_u_2 <- multiomicCCA::generate_sbm_orthogonal(B_mat, membership_vec, centered = T)
  
  p_1 <- 20; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_1)*c(0.75,0.5)
  svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, K-1)
  svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, K-1)
  
  dat <- multiomicCCA::generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2,
                       noise_val = 2)
  list(dat = dat, true_membership_vec = true_membership_vec)
}