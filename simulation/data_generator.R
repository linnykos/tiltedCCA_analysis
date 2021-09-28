simulation_all <- function(idx){
  if(idx == 1) { simulation1() 
  } else if(idx == 2){ simulation2() 
  } else if(idx == 3){ simulation3() 
  } else if(idx == 4){ simulation4() 
  } else if(idx == 5){ simulation5() 
  } else if(idx == 6){ simulation6() 
  } else if(idx == 7){ simulation7() 
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
  K <- ncol(B_mat); 
  rho1 <- 1; rho2 <- 0.3
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  n <- length(membership_vec); true_membership_vec <- membership_vec
  svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(rho1*B_mat, membership_vec, centered = T)
  svd_u_2 <- multiomicCCA::generate_sbm_orthogonal(rho2*B_mat, membership_vec, centered = T)
  
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

# 3 clusters
simulation4 <- function(){
  B_mat <- matrix(c(0.9, 0.9, 0.1,
                    0.9, 0.9, 0.1,
                    0.1, 0.1, 0.9), 3, 3, byrow = T)
  K <- ncol(B_mat); n_clust <- 100
  
  true_membership_vec <- rep(1:3, each = n_clust)
  n <- length(true_membership_vec)
  
  membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
  svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(B_mat, membership_vec, centered = T)
  membership_vec <- c(rep(1, n_clust), rep(3, n_clust), rep(2, n_clust))
  svd_u_2 <- multiomicCCA::generate_sbm_orthogonal(B_mat, membership_vec, centered = T)

  p_1 <- 40; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(5,1); svd_d_2 <- sqrt(n*p_2)*c(5,1)
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

# same setting as above, but now w/ only 4 variables
simulation6 <- function(){
  B_mat <- matrix(c(0.9, 0.9, 0.3,  
                    0.9, 0.9, 0.3,
                    0.3, 0.3, 0.9), 3, 3, byrow = T)
  K <- ncol(B_mat); n_clust <- 100
  true_membership_vec <- rep(1:4, each = n_clust)
  n <- length(true_membership_vec)
  
  membership_vec <-  c(rep(1, n_clust), rep(2, n_clust), rep(3, 2*n_clust))
  svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(B_mat, membership_vec, centered = T)
  membership_vec <-  c(rep(1, n_clust), rep(3, n_clust), rep(2, n_clust), rep(3, n_clust))
  svd_u_2 <- multiomicCCA::generate_sbm_orthogonal(B_mat, membership_vec, centered = T)
  
  p_1 <- 40; p_2 <- 40
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
  svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, K-1)
  svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, K-1)
  
  dat <- multiomicCCA::generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2,
                                     noise_val = 2)
  list(dat = dat, true_membership_vec = true_membership_vec)
}

# trajectory setting where the first modality has pseudotime and
# the second modality has branch information
simulation7 <- function(){
  n <- 400
  branch_vec <- c(rep(1, each = n/2), rep(2, each = n/4), rep(3, each = n/4))
  pseudotime_vec <- unlist(lapply(1:3, function(x){
    idx <- which(branch_vec == x)
    seq(0, 1, length.out = length(idx))
  }))
  
  tmp_mat_1 <- t(sapply(pseudotime_vec, function(x){
    c(rnorm(1, mean = x, sd = 0.05),
      rnorm(1, mean = runif(1, min = -x, max = x), sd = 0.05))
  }))
  centers <- matrix(c(0,0, 0.3,1, -0.3,1, 1,0), 
                    nrow = 4, ncol = 2, byrow = T)
  tmp_idx <- t(sapply(1:n, function(i){
    branch <- branch_vec[i]
    pseudotime <- pseudotime_vec[i]
    if(branch == 1){
      if(pseudotime <= 0.5){
        sample(c(1,4), size = 1, prob = c(0.8, 0.2))
      } else {
        sample(c(1,4), size = 1, prob = c(0.2, 0.8))
      }
    } else {
      if(pseudotime <= 0.5){
        sample(c(1:3), size = 1, prob = c(1/3, 1/3, 1/3))
      } else {
        if(branch == 2){
          sample(c(1,2), size = 1, prob = c(0.2, 0.8))
        } else {
          sample(c(1,3), size = 1, prob = c(0.2, 0.8))
        }
      }
    }
  }))
  tmp_mat_2 <- centers[tmp_idx,]
  tmp_mat_2 <- tmp_mat_2 + rnorm(prod(dim(tmp_mat_2)), mean = 0, sd = 0.05)

  tmp_mat_1 <- scale(tmp_mat_1, center = T, scale = T)
  tmp_mat_2 <- scale(tmp_mat_2, center = T, scale = T)
  
  svd_u_1 <- svd(tmp_mat_1)$u
  svd_u_2 <-svd(tmp_mat_2)$u
  
  p_1 <- 20; p_2 <- 20
  svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_1)*c(0.75,0.5)
  K <- 2
  svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, K)
  svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, K)
  
  dat <- multiomicCCA::generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2,
                                     noise_val = 2)
  list(dat = dat, true_branch = branch_vec,
       true_pseudotime = pseudotime_vec)
}