rm(list=ls())
library(simulator)
library(multiomicCCA)
source("../multiomicCCA_analysis/simulation/data_generator.R")

df_param <- data.frame(setting = c(1, 2, 3, 4, 5, 6), 
                       rank_1 =  c(2, 2, 2, 2, 2, 2),
                       rank_2 =  c(2, 2, 2, 2, 2, 2))

rule <- function(vec){
  simulation_all(vec$setting)
}

criterion <- function(dat, vec, y){
  dcca_res <- dcca_factor(dat$dat$mat_1, dat$dat$mat_2, rank_1 = vec$rank_1, rank_2 = vec$rank_2, 
                          apply_shrinkage = F, verbose = F)
  dcca_decomp <- dcca_decomposition(dcca_res, rank_c = min(vec$rank_1, vec$rank_2), verbose = F)
  
  list(dat = dat$dat, true_membership_vec = dat$true_membership_vec,
       dcca_decomp = dcca_decomp)
}

##################

res <- simulator::simulator(rule, criterion, df_param, ntrials = 1,
                            cores = 1, verbose = T)

###################

# choose a particular setting to investigate
i <- 1
dcca_decomp <- res[[i]][[1]]$result$dcca_decomp
true_membership_vec <- as.factor(res[[i]][[1]]$result$true_membership_vec)

dcca_decomp$distinct_perc_2
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec)
par(mar = c(4,4,4,0.5))
plot_data(dcca_decomp, membership_vec = true_membership_vec, pca = T)

par(mar = c(4,4,4,4), mfrow = c(1,1))
plot_summary(dcca_decomp)

par(mar = c(4,4,4,0.5))
plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec)
par(mar = c(4,4,4,0.5))
plot_scores_heatmap(dcca_decomp, membership_vec = true_membership_vec, scaling_power = 0.5)

par(mar = c(4,4,4,0.5))
plot_scores(dcca_decomp, membership_vec = true_membership_vec, decomposition = T)

par(mar = c(4,4,4,0.5), mfrow = c(1,2))
plot_decomposition_2d(dcca_decomp$score_1[,1], dcca_decomp$score_2[,1], dcca_decomp$common_score[,1], 
                      main = bquote(atop(bold("Canonical dim. 1"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[1],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$distinct_perc_2[1],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")
plot_decomposition_2d(dcca_decomp$score_1[,2], dcca_decomp$score_2[,2], dcca_decomp$common_score[,2], 
                      main = bquote(atop(bold("Canonical dim. 2"), ~ lambda ~ plain(":") ~ .(round(dcca_decomp$cca_obj[2],2)) 
                                         ~ plain(",")~ gamma ~ plain(":") ~ .(round(dcca_decomp$distinct_perc_2[2],2)))),
                      xlab = "Dataset 1", ylab = "Dataset 2")

par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = F, main_addition = "\nDataset 1",
                add_noise = T)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = F, data_2 = T, main_addition = "\nDataset 2",
                add_noise = T)
par(mar = c(4,4,4,0.5))
plot_embeddings(dcca_decomp, true_membership_vec, data_1 = T, data_2 = T, main_addition = "\nBoth datasets",
                add_noise = T)

set.seed(10)
clisi_1 <- clisi_information(dcca_decomp$common_mat_1, dcca_decomp$distinct_mat_1,
                             true_membership_vec, rank_c = 2, rank_d = 2, 
                             nn = round(0.5*min(table(true_membership_vec))),
                             radius_quantile = 0.9, max_subsample_clisi = min(table(true_membership_vec)),
                             verbose = F)
clisi_2 <- clisi_information(dcca_decomp$common_mat_2, dcca_decomp$distinct_mat_2,
                             true_membership_vec, rank_c = 2, rank_d = 2, 
                             nn = round(0.5*min(table(true_membership_vec))),
                             radius_quantile = 0.9, max_subsample_clisi = min(table(true_membership_vec)),
                             verbose = F)
plot_clisi(clisi_1, clisi_2)
plot_clisi_legend(clisi_1)


######

for(i in 1:6){
  dcca_decomp <- res[[i]][[1]]$result$dcca_decomp
  true_membership_vec <- as.factor(res[[i]][[1]]$result$true_membership_vec)
  set.seed(10)
  clisi_1 <- clisi_information(dcca_decomp$common_mat_1, dcca_decomp$distinct_mat_1,
                               true_membership_vec, rank_c = 2, rank_d = 2, 
                               nn = round(0.5*min(table(true_membership_vec))),
                               radius_quantile = 0.9, max_subsample_clisi = min(table(true_membership_vec)),
                               verbose = F)
  clisi_2 <- clisi_information(dcca_decomp$common_mat_2, dcca_decomp$distinct_mat_2,
                               true_membership_vec, rank_c = 2, rank_d = 2, 
                               nn = round(0.5*min(table(true_membership_vec))),
                               radius_quantile = 0.9, max_subsample_clisi = min(table(true_membership_vec)),
                               verbose = F)
  png(paste0("../../out/simulation/Writeup13/Writeup13_simulation", i, "_clisi.png"), height = 900, width = 1500, units = "px", res = 300)
  plot_clisi(clisi_1, clisi_2, col_vec = 1:length(levels(true_membership_vec)))
  graphics.off()
}

