rm(list=ls())
library(multiomicCCA)
source("../multiomicCCA_analysis/simulation/data_generator.R")
source("../multiomicCCA_analysis/simulation/plotting.R")
source("../multiomicCCA_analysis/simulation/wnn.R")
source("../multiomicCCA_analysis/simulation/pca_combine.R")
source("../multiomicCCA_analysis/simulation/plotting.R")
source("../multiomicCCA_analysis/simulation/dcca_custom.R")
source("../multiomicCCA_analysis/simulation/jive.R")

df_param <- data.frame(setting = c(1, 4, 3, 6, 7),
                       fix_tilt_perc = c(NA, NA, NA, NA, NA))

for(row_idx in 1:4){
  print(row_idx)
  set.seed(10)
  setting <- df_param[row_idx, "setting"]
  fix_tilt_perc <- df_param[row_idx, "fix_tilt_perc"]
  if(is.na(fix_tilt_perc)){
    fix_tilt_perc <- F
  }
  dat <- simulation_all(setting)
  
  set.seed(10)
  dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, 
                          dims_1 = 1:2, dims_2 = 1:2, 
                          fix_tilt_perc = fix_tilt_perc,
                          center_1 = F, center_2 = F,
                          scale_1 = F, scale_2 = F,
                          verbose = F,
                          metacell_clustering = as.factor(dat$true_membership_vec))
  print(dcca_res$tilt_perc)
  print(dcca_res$df_percentage)
  svd_1 <- dcca_res$svd_1
  svd_2 <- dcca_res$svd_2
  dims_1 <- 1:2
  dims_2 <- 1:2
  mat_1 <- dat$mat_1
  mat_2 <- dat$mat_2
  dimred1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)
  dimred2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d)
  dimred_combined <- pca_combine(svd_1, svd_2, dims_1, dims_2)
  dimred_wnn <- wnn_custom(mat_1, mat_2, svd_1, svd_2, dims_1, dims_2)
  dimred_dcca <- dcca_custom(dcca_res)
  tmp <- jive(mat_1, mat_2, r = 2)
  dimred_jive <- tmp$embedding
  
  if(setting < 7){
    color_vec <- scales::hue_pal()(length(unique(dat$true_membership_vec)))
    
    # dataset 1
    plot_custom_cluster(dimred1,
                        membership_vec = dat$true_membership_vec,
                        color_vec = color_vec,
                        filename = paste0("../../out/figures/Writeup14d_simulation/simulation", setting,
                                          "_embedding1.png"),
                        main = "Dataset 1 (PCA)")
    
    # dataset 2
    plot_custom_cluster(dimred2,
                        membership_vec = dat$true_membership_vec,
                        color_vec = color_vec,
                        filename = paste0("../../out/figures/Writeup14d_simulation/simulation", setting,
                                          "_embedding2.png"),
                        main = "Dataset 2 (PCA)")
    
    # pca combined
    plot_custom_cluster(dimred_combined,
                        membership_vec = dat$true_membership_vec,
                        color_vec = color_vec,
                        filename = paste0("../../out/figures/Writeup14d_simulation/simulation", setting,
                                          "_pcacombined.png"),
                        main = "PCA Combined")
    
    # wnn
    plot_custom_cluster(dimred_wnn,
                        membership_vec = dat$true_membership_vec,
                        color_vec = color_vec,
                        filename = paste0("../../out/figures/Writeup14d_simulation/simulation", setting,
                                          "_wnn.png"),
                        main = "WNN")
    
    # jive
    plot_custom_cluster(dimred_jive,
                        membership_vec = dat$true_membership_vec,
                        color_vec = color_vec,
                        filename = paste0("../../out/figures/Writeup14d_simulation/simulation", setting,
                                          "_jive.png"),
                        main = "JIVE")
    
    xlim <- range(c(dimred_dcca$common[,1], dimred_dcca$distinct_1[,1], dimred_dcca$distinct_2[,1]))
    ylim <- range(c(dimred_dcca$common[,2], dimred_dcca$distinct_1[,2], dimred_dcca$distinct_2[,2]))
    
    plot_custom_cluster(dimred_dcca$common,
                        membership_vec = dat$true_membership_vec,
                        color_vec = color_vec,
                        filename = paste0("../../out/figures/Writeup14d_simulation/simulation", setting,
                                          "_dcca_common.png"),
                        main = paste0("Tilted-CCA (Common): Tilt of ", dcca_res$tilt_perc),
                        xlim = xlim, ylim = ylim)
    
    plot_custom_cluster(dimred_dcca$distinct_1,
                        membership_vec = dat$true_membership_vec,
                        color_vec = color_vec,
                        filename = paste0("../../out/figures/Writeup14d_simulation/simulation", setting,
                                          "_dcca_distinct_1.png"),
                        main = "Tilted-CCA (Distinct 1)",
                        xlim = xlim, ylim = ylim)
    
    plot_custom_cluster(dimred_dcca$distinct_2,
                        membership_vec = dat$true_membership_vec,
                        color_vec = color_vec,
                        filename = paste0("../../out/figures/Writeup14d_simulation/simulation", setting,
                                          "_dcca_distinct_2.png"),
                        main = "Tilted-CCA (Distinct 2)",
                        xlim = xlim, ylim = ylim)
    
  } else {
    
  }
}

