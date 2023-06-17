###########

k <- 11
mat_1 <- result_list[[k]]$mat_1
mat_2 <- result_list[[k]]$mat_2
true_cluster <- result_list[[k]]$true_cluster
par(mfrow = c(1,2))
tmp <- svd_func(mat_1)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "Data 1",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
tmp <- svd_func(mat_2)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "Data 2",
     xlab = "", ylab = "",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)

###########

k <- 11
cex <- 1.3; cex.main <- 1.5
multiSVD_obj <- result_list[[k]]$example_res
true_cluster <- result_list[[k]]$true_cluster
par(mfrow = c(1,3))
tmp <- svd_func(multiSVD_obj$tcca_obj$common_score)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "Common embedding",
     xlab = "Common's dim. 1", ylab = "Common's dim. 2",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T,
     cex.lab = cex, cex = cex, cex.axis = cex, cex.main = cex.main)
tmp <- svd_func(multiSVD_obj$tcca_obj$distinct_score_1)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "Modality 1's distinct embedding",
     xlab = "Distinct-1's dim. 1", ylab = "Distinct-1's dim. 2",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T,
     cex.lab = cex, cex = cex, cex.axis = cex, cex.main = cex.main)
tmp <- svd_func(multiSVD_obj$tcca_obj$distinct_score_2)
plot(tmp[plot_idx,1], tmp[plot_idx,2],
     main = "Modality 2's distinct embedding",
     xlab = "Distinct-2's dim. 1", ylab = "Distinct-2's dim. 2",
     pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T,
     cex.lab = cex, cex = cex, cex.axis = cex, cex.main = cex.main)
