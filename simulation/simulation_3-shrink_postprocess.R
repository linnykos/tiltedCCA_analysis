rm(list=ls())
load("../../out/simulation/simulation_3-shrink.RData")

alignment_mat <- matrix(NA, nrow = 3, ncol = length(result_list))
colnames(alignment_mat) <- names(result_list)
rownames(alignment_mat) <- paste0("quantile_", c(25,50,75))
for(j in 1:length(result_list)){
  alignment_mat[,j] <- stats::quantile(
    result_list[[j]]$alignment,
    probs = c(0.25,0.5,0.75)
  )
}

table_mat <- matrix(NA, nrow = 3, ncol = length(result_list))
colnames(table_mat) <- names(result_list)
rownames(table_mat) <- 1:3
for(j in 1:length(result_list)){
  table_mat[,j] <- sapply(1:3, function(k){
    length(which(result_list[[j]]$clusters == k))
  })
  table_mat[,j] <- 100*table_mat[,j]/sum(table_mat[,j])
}


############################

# Grouped barplot
# see https://r-graph-gallery.com/211-basic-grouped-or-stacked-barplot.html
# and https://www.statology.org/r-height-must-be-vector-or-matrix/
png("../tiltedCCA_analysis/simulation/fig/simulation_3-shrink_clusters.png", 
    width = 2000, height = 1000, units = "px", res = 300)
barplot(table_mat,
        col=colors()[c(23,89,12)] , 
        border="white", 
        space=0.04,
        xlab="",
        axisnames = F)
graphics.off()

# alignment plot
png("../tiltedCCA_analysis/simulation/fig/simulation_3-shrink_alignment.png", 
    width = 2000, height = 1000, units = "px", res = 300)
ylim <- range(alignment_mat)
x_vec <- seq(0, 1, length.out=11)
plot(NA, xlim = range(x_vec), ylim = ylim, 
     xlab = "", ylab = "", 
     main = "", bty = "n")
axis(1)
axis(2)
polygon(x = c(x_vec, rev(x_vec)),
        y = c(alignment_mat[1,], rev(alignment_mat[3,])),
        col = rgb(0.5, 0.5, 0.5, 0.5),
        border = "black")
points(x = x_vec, y = alignment_mat[2,], pch = 16, cex = 2)
lines(x = x_vec, y = alignment_mat[2,], lwd = 2, lty = 2)
graphics.off()

############################
source("../tiltedCCA_analysis/simulation/simulation_3-shrink_functions.R")

set.seed(10)
plot_idx <- sample(1:300) # 300 is equal to nrow(mat_1)
orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
blue_col <- rgb(129, 139, 191, maxColorValue = 255)
col_vec <- c(purple_col, orange_col, blue_col)

for(k in c(1,6,11)){
  mat_2 <- result_list[[k]]$mat_2
  tmp <- svd_func(mat_2)
  
  true_cluster <- result_list[[k]]$true_cluster
  png(paste0("../tiltedCCA_analysis/simulation/fig/simulation_3-mat-2_", names(result_list)[k], ".png"), 
      width = 800, height = 800, units = "px", res = 300)
  par(mar = c(3,3,0.5,0.5))
  plot(tmp[plot_idx,1], tmp[plot_idx,2],
       main = "",
       xlab = "", ylab = "",
       pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
  graphics.off()
}

###########

for(k in c(1,6,11)){
  multiSVD_obj <- result_list[[k]]$example_res
  true_cluster <- result_list[[k]]$true_cluster
  tmp <- svd_func(multiSVD_obj$tcca_obj$common_score)
  
  png(paste0("../tiltedCCA_analysis/simulation/fig/simulation_3-common-embedding_", names(result_list)[k], ".png"), 
      width = 800, height = 800, units = "px", res = 300)
  par(mar = c(3,3,0.5,0.5))
  plot(tmp[plot_idx,1], tmp[plot_idx,2],
       main = "",
       xlab = "", ylab = "",
       pch = 16, col = col_vec[true_cluster[plot_idx]], asp = T)
  graphics.off()
}
