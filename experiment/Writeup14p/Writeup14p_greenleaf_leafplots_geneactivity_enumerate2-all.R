rm(list=ls())
library(Seurat); library(Signac)
load("../../../../out/main/10x_greenleaf_tcca_RNA-geneActivity.RData")

var_features <- colnames(multiSVD_obj$common_mat_1)
var_features <- var_features[which(paste0("ATAC-",var_features) %in% colnames(multiSVD_obj$common_mat_2))]
var_features <- sort(var_features)

cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_mat <- multiSVD_obj$common_mat_1[cell_idx[order_vec],var_features]
atac_mat <- multiSVD_obj$common_mat_2[cell_idx[order_vec],paste0("ATAC-", var_features)]

n <- nrow(rna_mat)
p <- ncol(rna_mat)
pred_rna_mat <- sapply(1:p, function(j){
  print(j)
  
  tmp_df <- data.frame(y = rna_mat[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
pred_atac_mat <- sapply(1:p, function(j){
  print(j)
  
  tmp_df <- data.frame(y = atac_mat[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
colnames(pred_atac_mat) <- colnames(atac_mat)
colnames(pred_rna_mat) <- colnames(rna_mat)

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
color_vec <- rev(grDevices::colorRampPalette(base_palette)(n))

for(k in 1:ceiling(ncol(rna_mat)/30)){
  print(k)
  genes <- colnames(rna_mat)[((k-1)*30+1):min((k*30), ncol(rna_mat))]
  
  png(paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_greenleaf_leafplots_RNA-geneActivity_common_enumerate_", 
             k, ".png"),
      height = 3500, width = 3000, units = "px", res = 300)
  par(mfrow = c(6,5), mar = c(4,4,4,0.5))
  for(gene in genes){
    set.seed(10)
    idx <- sample(1:n)
    
    x_vec1 <- atac_mat[,paste0("ATAC-", gene)]
    x_vec2 <- pred_atac_mat[,paste0("ATAC-", gene)]
    y_vec1 <- rna_mat[,gene]
    y_vec2 <- pred_rna_mat[,gene]
    
    x_mean <- mean(x_vec2); x_vec1 <- x_vec1- x_mean; x_vec2 <- x_vec2 - x_mean
    x_sd <- sd(x_vec2); x_vec1 <- x_vec1/x_sd; x_vec2 <- x_vec2/x_sd
    y_mean <- mean(y_vec2); y_vec1 <- y_vec1- y_mean; y_vec2 <- y_vec2 - y_mean
    y_sd <- sd(y_vec2); y_vec1 <- y_vec1/y_sd; y_vec2 <- y_vec2/y_sd
    
    plot(x = x_vec1[idx], y = y_vec1[idx],
         pch = 16, col = color_vec[idx], main = gene, cex = 0.5,
         xlim = quantile(x_vec1, probs = c(0.01, 0.99)),
         ylim = quantile(y_vec1, probs = c(0.01, 0.99)),
         xlab = "ATAC", ylab = "RNA")
    
    points(x = x_vec2[idx], 
           y = y_vec2[idx], 
           pch = 16, col = "white", cex = 2)
    points(x = x_vec2[idx], 
           y = y_vec2[idx], 
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
    lines(c(x_right,x_left), c(y_top, y_bot), lwd = 2, lty = 2)
  }
  graphics.off()
}