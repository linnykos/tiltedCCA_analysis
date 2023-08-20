rm(list=ls())
library(tiltedCCA)
library(mclust)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

source("../tiltedCCA_analysis/simulation/simulation_3-shrink_functions.R")

trial <- 1

# simulate data
for(shrink in c(0, 0.5, 1)){
  print(paste0("Working on shrink ", shrink))
  
  set.seed(10*trial)
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
                                          latent_k = 10,
                                          num_neigh = 3,
                                          bool_cosine = T,
                                          bool_intersect = F,
                                          min_deg = 0)
  multiSVD_obj <- tiltedCCA::tiltedCCA(input_obj = multiSVD_obj,
                                       fix_tilt_perc = F)
  multiSVD_obj <- tiltedCCA::fine_tuning(input_obj = multiSVD_obj,
                                         verbose = 0)
  
  ###########
  
  downsample_idx <- c(1:40, 101:140, 201:240)
  png(paste0("../../out/figures/simulation/simulation_3_target-manifold_shrink-", shrink, ".png"), 
      width = 400, height = 400, units = "px", res = 300)
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  image(as.matrix(multiSVD_obj$snn_list$common_snn)[downsample_idx,downsample_idx], asp = T,
        col = c(rgb(230, 230, 230, maxColorValue = 255),
                rgb(223, 83, 107, maxColorValue = 255)),
        xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  for(i in 1:2){
    lines(x = c(0,1), y = rep(i/3, 2), lty = 2)
    lines(x = rep(i/3, 2), y = c(0,1), lty = 2)
  }
  graphics.off()
  
  commm_snn <- tiltedCCA:::.form_snn_mat(bool_cosine = multiSVD_obj$param$snn_bool_cosine,
                                         bool_intersect = multiSVD_obj$param$snn_bool_intersect,
                                         mat = multiSVD_obj$tcca_obj$common_score, 
                                         min_deg = multiSVD_obj$param$snn_min_deg,
                                         num_neigh = multiSVD_obj$param$snn_num_neigh,
                                         verbose = 0)
  
  png(paste0("../../out/figures/simulation/simulation_3_common-graph_shrink-", shrink, ".png"), 
      width = 400, height = 400, units = "px", res = 300)
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  image(as.matrix(commm_snn)[downsample_idx,downsample_idx], asp = T,
        col = c(rgb(230, 230, 230, maxColorValue = 255),
                rgb(223, 83, 107, maxColorValue = 255)),
        xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  for(i in 1:2){
    lines(x = c(0,1), y = rep(i/3, 2), lty = 2)
    lines(x = rep(i/3, 2), y = c(0,1), lty = 2)
  }
  graphics.off()
}
