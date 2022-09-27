rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_greenleaf_tcca_RNA-geneActivity.RData")

gene_ordering <- read.csv("../../../out/main/10x_greenleaf_developmentalGenes_heatmapOrder_gene_unformatted.txt",
                          header = F)

cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)

rna_common <- multiSVD_obj$common_mat_1[cell_idx[order_vec],gene_ordering[,1]]
atac_common <- multiSVD_obj$common_mat_2[cell_idx[order_vec],paste0("ATAC-", gene_ordering[,1])]

construction_smoothed_matrix <- function(mat){
  n <- nrow(mat)
  p <- ncol(mat)
  
  pred_mat <- sapply(1:p, function(j){
    print(j)
    
    tmp_df <- data.frame(y = mat[,j], x = 1:n)
    reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
    x_vec <- 1:n
    y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
    y_vec$Estimation[,"Pred"]
  })
  
  colnames(pred_mat) <- colnames(mat)
  pred_mat
}

pred_rna_common <- construction_smoothed_matrix(rna_common)
pred_atac_common <- construction_smoothed_matrix(atac_common)

n <- nrow(pred_rna_common)
color_vec <- rev(grDevices::colorRampPalette(c(rgb(191, 74, 223, maxColorValue = 255),
                                               rgb(0.8, 0.8, 0.8),
                                               rgb(239, 158, 88, maxColorValue = 255)))(n))
example_color <- rgb(1,1,1,0.5)
color_vec_trans <- sapply(color_vec, function(x){
  paste0(x, substr(example_color, start = 8, stop = 10))
})

for(k in 1:2){
  png(paste0("../../../out/figures/main/10x_greenleaf_leafplot_mainplot_seuratGActivity_common_appendix_all_", k, ".png"),
      height = 4000, width = 4000, units = "px", res = 500)
  par(mar = c(2.5,2.5,3,0.5), mfrow = c(5,5))
  
  for(gene in colnames(rna_common)[((k-1)*25+1):(25*k)]){
    set.seed(10)
    idx <- sample(1:n)
    
    x_vec1 <- atac_common[,paste0("ATAC-", gene)]
    x_vec2 <- pred_atac_common[,paste0("ATAC-", gene)]
    y_vec1 <- rna_common[,gene]
    y_vec2 <- pred_rna_common[,gene]
    
    x_mean <- mean(x_vec2); x_vec1 <- x_vec1- x_mean; x_vec2 <- x_vec2 - x_mean
    x_sd <- sd(x_vec2); x_vec1 <- x_vec1/x_sd; x_vec2 <- x_vec2/x_sd
    y_mean <- mean(y_vec2); y_vec1 <- y_vec1- y_mean; y_vec2 <- y_vec2 - y_mean
    y_sd <- sd(y_vec2); y_vec1 <- y_vec1/y_sd; y_vec2 <- y_vec2/y_sd
    
    plot(x = x_vec1[idx], y = y_vec1[idx],
         pch = 16, col = color_vec_trans[idx], cex = 0.75,
         xlim = quantile(x_vec1, probs = c(0.005, 0.995)),
         ylim = quantile(y_vec1, probs = c(0.005, 0.995)),
         xlab = "", ylab = "",
         xaxt = "n", yaxt = "n", bty = "n",
         main = gene)
    axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
    axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
    
    points(x = x_vec2[idx], 
           y = y_vec2[idx], 
           pch = 16, col = "white", cex = 2)
    points(x = rev(x_vec2), 
           y = rev(y_vec2), 
           pch = 16, col = rev(color_vec), cex = 1)
    
    points(x = x_vec2[1], 
           y = y_vec2[1],
           pch = 16, col = "white", cex = 4)
    points(x = x_vec2[1], 
           y = y_vec2[1], 
           pch = 16, col = color_vec[1], cex = 3)
    
    x_right <- mean(x_vec1[x_vec1 >= quantile(x_vec1, probs = 0.9)])
    x_left <- mean(x_vec1[x_vec1 <= quantile(x_vec1, probs = 0.1)])
    y_top <- mean(y_vec1[y_vec1 >= quantile(y_vec1, probs = 0.9)])
    y_bot <- mean(y_vec1[y_vec1 <= quantile(y_vec1, probs = 0.1)])
    
    slope <- (y_top-y_bot)/(x_right-x_left)
    intercept <- y_top - slope*x_right
    x_left2 <- -1e4; x_right2 <- 1e4
    y_bot2 <- slope*x_left2+intercept; y_top2 <- slope*x_right2+intercept
    lines(c(x_right2,x_left2), c(y_top2, y_bot2), lwd = 2, lty = 2)
  }
  
  graphics.off()
}

####################

gene_vec <- c("DACH1", "FIGN", "NEUROD2", "CUX2", "SATB2", "SNTG2", 
              "GABRB2", "DOK5", "NYAP2", "SLIT1",  "GNAL", "CHL1", "SYT16",
              "ELMOD1", "SYNDIG1")

png(paste0("../../../out/figures/main/10x_greenleaf_leafplot_mainplot_seuratGActivity_common_appendix.png"),
    height = 4000, width = 2500, units = "px", res = 500)
par(mar = c(2.5,2.5,3,0.5), mfrow = c(5,3))

for(gene in gene_vec){
  set.seed(10)
  idx <- sample(1:n)
  
  x_vec1 <- atac_common[,paste0("ATAC-", gene)]
  x_vec2 <- pred_atac_common[,paste0("ATAC-", gene)]
  y_vec1 <- rna_common[,gene]
  y_vec2 <- pred_rna_common[,gene]
  
  x_mean <- mean(x_vec2); x_vec1 <- x_vec1- x_mean; x_vec2 <- x_vec2 - x_mean
  x_sd <- sd(x_vec2); x_vec1 <- x_vec1/x_sd; x_vec2 <- x_vec2/x_sd
  y_mean <- mean(y_vec2); y_vec1 <- y_vec1- y_mean; y_vec2 <- y_vec2 - y_mean
  y_sd <- sd(y_vec2); y_vec1 <- y_vec1/y_sd; y_vec2 <- y_vec2/y_sd
  
  plot(x = x_vec1[idx], y = y_vec1[idx],
       pch = 16, col = color_vec_trans[idx], cex = 0.75,
       xlim = quantile(x_vec1, probs = c(0.005, 0.995)),
       ylim = quantile(y_vec1, probs = c(0.005, 0.995)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", bty = "n",
       main = gene)
  axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
  axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
  
  points(x = x_vec2[idx], 
         y = y_vec2[idx], 
         pch = 16, col = "white", cex = 2)
  points(x = rev(x_vec2), 
         y = rev(y_vec2), 
         pch = 16, col = rev(color_vec), cex = 1)
  
  points(x = x_vec2[1], 
         y = y_vec2[1],
         pch = 16, col = "white", cex = 4)
  points(x = x_vec2[1], 
         y = y_vec2[1], 
         pch = 16, col = color_vec[1], cex = 3)
  
  x_right <- mean(x_vec1[x_vec1 >= quantile(x_vec1, probs = 0.9)])
  x_left <- mean(x_vec1[x_vec1 <= quantile(x_vec1, probs = 0.1)])
  y_top <- mean(y_vec1[y_vec1 >= quantile(y_vec1, probs = 0.9)])
  y_bot <- mean(y_vec1[y_vec1 <= quantile(y_vec1, probs = 0.1)])
  
  slope <- (y_top-y_bot)/(x_right-x_left)
  intercept <- y_top - slope*x_right
  x_left2 <- -1e4; x_right2 <- 1e4
  y_bot2 <- slope*x_left2+intercept; y_top2 <- slope*x_right2+intercept
  lines(c(x_right2,x_left2), c(y_top2, y_bot2), lwd = 2, lty = 2)
}

graphics.off()
