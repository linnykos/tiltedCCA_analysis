rm(list=ls())
library(tiltedCCA)
library(mclust)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

source("../tiltedCCA_analysis/simulation/simulation_1-dampeningB_functions.R")

################################
# Step 1: Generate the data, where Modality 1 is informative, but not Modality 2
################################

plot_idx <- sample(1:300) # 300 is equal to nrow(mat_1)
orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
blue_col <- rgb(129, 139, 191, maxColorValue = 255)
col_vec <- c(purple_col, orange_col, blue_col)

################################
# Step 2: Apply Tilted-CCA for sequence of separations: 100 trials for 5 different values
################################

trials <- 50

shrink_vec <- seq(0, 1, length.out=11)
len <- length(shrink_vec)
result_list <- lapply(1:len, function(i){
  list(clustering_1 = NA,
       clustering_2 = NA,
       mat_1 = NA,
       mat_2 = NA,
       example_res = NA, 
       true_cluster = NA,
       clusters = rep(NA, trials),
       alignment = rep(NA, trials))
})
names(result_list) <- paste0("shrink_", shrink_vec)

for(shrink in shrink_vec){
  print("=====")
  print(paste0("Working on shrink ", round(shrink,2)))
  for(trial in 1:trials){
    set.seed(10*trial)
    if(trial %% floor(trials/10) == 0) cat('*')
    
    # simulate data
    res <- simulation_function(shrink = shrink)
    
    set.seed(10)
    multiSVD_obj <- tiltedCCA::create_multiSVD(mat_1 = res$mat_1, mat_2 = res$mat_2,
                                               dims_1 = 1:3, dims_2 = 1:3,
                                               center_1 = T, center_2 = T,
                                               normalize_row = T,
                                               normalize_singular_value = F,
                                               recenter_1 = F, recenter_2 = F,
                                               rescale_1 = F, rescale_2 = F,
                                               scale_1 = T, scale_2 = T)
    multiSVD_obj <- tiltedCCA::form_metacells(input_obj = multiSVD_obj,
                                              large_clustering_1 = res$clustering_1, 
                                              large_clustering_2 = res$clustering_2,
                                              num_metacells = NULL)
    multiSVD_obj <- tiltedCCA::compute_snns(input_obj = multiSVD_obj,
                                            latent_k = 3,
                                            num_neigh = 20,
                                            bool_cosine = T,
                                            bool_intersect = F,
                                            min_deg = 0)
    multiSVD_obj <- tiltedCCA::tiltedCCA(input_obj = multiSVD_obj,
                                         fix_tilt_perc = F)
    multiSVD_obj <- tiltedCCA::fine_tuning(input_obj = multiSVD_obj,
                                           verbose = 0)
    multiSVD_obj <- tiltedCCA::tiltedCCA_decomposition(multiSVD_obj)
    
    if(trial == 1){
      result_list[[paste0("shrink_", shrink)]]$clustering_1 <- res$clustering_1
      result_list[[paste0("shrink_", shrink)]]$clustering_2 <- res$clustering_2
      result_list[[paste0("shrink_", shrink)]]$mat_1 <- res$mat_1
      result_list[[paste0("shrink_", shrink)]]$mat_2 <- res$mat_2
      result_list[[paste0("shrink_", shrink)]]$true_cluster <- res$true_cluster
      result_list[[paste0("shrink_", shrink)]]$example_res <- multiSVD_obj
    }
    
    tmp <- svd_func(multiSVD_obj$tcca_obj$common_score)
    
    mclust_res <- mclust::Mclust(tmp, 
                                 G = 1:3,
                                 modelNames = "VVV",
                                 verbose = F)
    result_list[[paste0("shrink_", shrink)]]$clusters[trial] <- mclust_res$G
    result_list[[paste0("shrink_", shrink)]]$alignment[trial] <- compute_alignment(
      mat_1 = res$mat_1,
      mat_2 = res$mat_2,
      multiSVD_obj = multiSVD_obj,
      true_cluster = res$true_cluster
    )
  }
  
  print("\n")
  print(table(result_list[[paste0("shrink_", shrink)]]$clusters))
}

save(result_list, date_of_run, session_info,
     file = "../../out/simulation/simulation_1-dampeningB.RData")

###########

for(i in 1:length(result_list)){
  print(table(result_list[[i]]$clusters))
}

for(i in 1:length(result_list)){
  print(median(result_list[[i]]$alignment))
}

###########################################

k <- 11



