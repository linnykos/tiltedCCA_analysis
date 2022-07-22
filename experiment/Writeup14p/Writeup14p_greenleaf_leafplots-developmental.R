rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/main/10x_greenleaf_timeseries.RData")
load("../../../../out/main/10x_greenleaf_developmentalGenes.RData")

colnames(rna_predicted_mat) <- colnames(rna_common)
colnames(atac_predicted_mat) <- colnames(atac_pred)

n <- nrow(rna_common)
base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
color_vec <- rev(grDevices::colorRampPalette(base_palette)(n))

for(k in 1:2){
  genes <- selection_res$selected_variables[((k-1)*25+1):min((k*25), length(selection_res$selected_variables))]
  
  png(paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_greenleaf_leafplots-developmental_",
             k, ".png"),
      height = 3000, width = 3000, units = "px", res = 300)
  par(mfrow = c(5,5), mar = c(4,4,4,0.5))
  for(gene in genes){
    set.seed(10)
    idx <- sample(1:n)
    
    x_vec1 <- atac_pred[idx,gene]
    x_vec2 <- atac_predicted_mat[idx,gene]
    y_vec1 <- rna_common[idx,gene]
    y_vec2 <- rna_predicted_mat[idx,gene]
    
    x_mean <- mean(x_vec2); x_vec1 <- x_vec1- x_mean; x_vec2 <- x_vec2 - x_mean
    x_sd <- sd(x_vec2); x_vec1 <- x_vec1/x_sd; x_vec2 <- x_vec2/x_sd
    y_mean <- mean(y_vec2); y_vec1 <- y_vec1- y_mean; y_vec2 <- y_vec2 - y_mean
    y_sd <- sd(y_vec2); y_vec1 <- y_vec1/y_sd; y_vec2 <- y_vec2/y_sd
    
    plot(x = x_vec1, y = y_vec1,
         pch = 16, col = color_vec[idx], main = gene, cex = 0.5,
         xlab = "ATAC", ylab = "RNA", asp = T)
    
    points(x = x_vec2, 
           y = y_vec2, 
           pch = 16, col = "white", cex = 2)
    points(x = x_vec2, 
           y = y_vec2, 
           pch = 16, col = color_vec[idx], cex = 1.5)
    
    points(x = x_vec2[1], 
           y = y_vec2[1], 
           pch = 16, col = "white", cex = 3)
    points(x = x_vec2[1], 
           y = y_vec2[1], 
           pch = 16, col = color_vec[1], cex = 2.5)
    
    x_right <- mean(x_vec1[x_vec1 >= quantile(x_vec1, probs = 0.9)])
    x_left <- mean(x_vec1[x_vec1 <= quantile(x_vec1, probs = 0.1)])
    y_top <- mean(y_vec1[y_vec1 >= quantile(y_vec1, probs = 0.9)])
    y_bot <- mean(y_vec1[y_vec1 <= quantile(y_vec1, probs = 0.1)])
    lines(c(2*x_right,2*x_left), c(2*y_top, 2*y_bot), lwd = 2, lty = 2)
  }
  graphics.off()
}