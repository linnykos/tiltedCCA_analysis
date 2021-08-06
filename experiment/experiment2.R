rm(list=ls())
library(customSimulator)
library(networkSoSD)

session_info <- devtools::session_info()
date_of_run <- Sys.time()
source_code_info <- readLines("../simulation/simulation_rho.R")
run_suffix <- ""

# df_param <- cbind(3, 500, 100, seq(0.025, 0.2, length.out = 15), 0.4, 0.1, 0.5)
df_param <- cbind(3, 30, 10, seq(0.025, 0.2, length.out = 2), 0.4, 0.1, 0.5)
colnames(df_param) <- c("K", "n", "L", "rho", "mem_prop1", "mem_prop2", "mem_prop3")
df_param <- as.data.frame(df_param)

# ntrials <- 50
ntrials <- 5
ncores <- 2
worker_variables <- list(lmfo = lmfo, coreg = coreg, vec = 1:10)

###############################

generator <- function(vec, worker_variables){
  .l2norm <- function(x){sqrt(sum(x^2))}
  vec1 <- c(1,1,sqrt(2))
  vec2 <- c(1,1,-sqrt(2))
  vec3 <- c(-1,1,0)
  eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
  B1 <- eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat)
  B2 <- eigen_mat %*% diag(c(1.5, 0.2, -0.4)) %*% t(eigen_mat)
  
  n <- as.numeric(vec["n"]); L <- as.numeric(vec["L"]); rho <- as.numeric(vec["rho"])
  K <- as.numeric(vec["K"])
  mem_prop1 <- as.numeric(vec["mem_prop1"]); mem_prop2 <- as.numeric(vec["mem_prop2"])
  mem_prop3 <- as.numeric(vec["mem_prop3"])
  
  membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
  if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))
  
  prob_mat1 <- networkSoSD::compute_prob_mat(rho*B1, membership_vec)
  prob_mat2 <- networkSoSD::compute_prob_mat(rho*B2, membership_vec)
  prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})
  
  adj_list <- lapply(1:L, function(i){networkSoSD::generate_adjaceny_mat(prob_list[[i]])})
  
  list(adj_list = adj_list)
}

executor <- function(dat, vec, y, worker_variables){
  K <- as.numeric(vec["K"])
  
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias")
  res1 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of adding
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "sum")
  res2 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of ss
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
  res3 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of flattening
  set.seed(10)
  flat_mat <- networkSoSD::flatten(dat$adj_list)
  res4 <- networkSoSD::spectral_clustering(flat_mat, K = K, weighted = F)
  
  ### now all the weighted versions
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias")
  res1b <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)
  
  # try naive method of adding
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "sum")
  res2b <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)
  
  # try naive method of ss
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
  res3b <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)
  
  # try naive method of flattening
  set.seed(10)
  flat_mat <- networkSoSD::flatten(dat$adj_list)
  res4b <- networkSoSD::spectral_clustering(flat_mat, K = K, weighted = T)
  
  ### now the greedy method
  set.seed(10)
  res5 <- networkSoSD::greedy_refinement(dat$adj_list, K = K)$cluster
  
  # yuguo's methods
  set.seed(10)
  res6 <- networkSoSD::lmfo(dat$adj_list, k = K)
  
  # set.seed(10)
  reg_val <- max(sapply(dat$adj_list, function(x){4*max(abs(RSpectra::eigs_sym(x, k = min(K,5))$values))}))
  tmp <- networkSoSD::coreg(dat$adj_list, k = K, beta = reg_val,
                            verbose = F, max_iter = 50)
  res7 <- tmp[[length(dat$adj_list)+1]]
  
  list(res_ss_debias_F = res1, res_ss_debias_T = res1b, 
       res_sum_F = res2, res_sum_T = res2b, 
       res_ss_F = res3, res_ss_T = res3b,
       res_flat_F = res4, res_flat_T = res4b,
       res_greedy = res5, chen_linked = res6, chen_coreg = res7)
}

i <- 1; y <- 1; set.seed(y); zz <- executor(generator(df_param[i,], worker_variables), df_param[i,], y, worker_variables); zz
