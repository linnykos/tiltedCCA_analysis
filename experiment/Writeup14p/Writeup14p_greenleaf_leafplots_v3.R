rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/Writeup14p/Writeup14p_10x_greenleaf_tcca.RData")

rna_mat <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_mat <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)

cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_mat <- rna_mat[cell_idx[order_vec],]
atac_mat <- atac_mat[cell_idx[order_vec],]

mentioned_genes <- c("ALDH2", "APOE", "AQP4", "ASCL1", "BHLHE22",
                     "BHLHE40", "C16orf89", "CAV2", "DLX2", "DOK5",
                     "DUSP1", "EOMES", "ETV4", "FOS", "FOXJ1", "GLI3",
                     "HAS2", "HES1", "HES4", "HSPA1A", "HSPA1B",
                     "ID3", "IGFBP7", "JUN", "KIF1A", "LIMCH1",
                     "MBP", "MEF2C", "NEUROD1", "NEUROD2", "NEUROD4",
                     "NEUROD6", "NEUROG1", "NEUROG2", "NFIA", "NFIB", "NFIC",
                     "NHLH1", "NR2F1", "PAX6", "RFX4", "RUNX1",
                     "OLIG1", "OLIG2", "SOX2", "SOX3", "SOX6",
                     "SOX9", "SOX10", "SOX21", "SPARCL1", "SNCB", "TBX",
                     "TNC", "TOP2A", "TRB1", "WNT11")
mentioned_genes <- setdiff(mentioned_genes, c("NFIA", "RUNX1", "SOX6"))
mentioned_genes <- sort(intersect(mentioned_genes, colnames(rna_mat)))

#####################################

rna_mat2 <- t(scale(t(rna_mat)))
atac_mat2 <- t(scale(t(atac_mat)))

diff_mat <- (rna_mat2 - atac_mat2)/abs(rna_mat2)
rna_mat2 <- rna_mat2[,mentioned_genes]
atac_mat2 <- atac_mat2[,mentioned_genes]
diff_mat <- diff_mat[,mentioned_genes]
diff_mat <- pmax(pmin(diff_mat, 15), -15)

p <- ncol(diff_mat); n <- nrow(diff_mat)
pred_rna_mat <- sapply(1:p, function(j){
  print(j)
  
  tmp_df <- data.frame(y = rna_mat2[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
colnames(pred_rna_mat) <- colnames(rna_mat2)
pred_diff_mat <- sapply(1:p, function(j){
  print(j)
  
  tmp_df <- data.frame(y = diff_mat[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
colnames(pred_diff_mat) <- colnames(diff_mat)

round(pred_diff_mat[round(seq(1,n,length.out=100)),"EOMES"], 2)
atac_mat3 <- pred_rna_mat
for(j in 1:p){
  atac_mat3[,j] <- pred_rna_mat[,j] - abs(pred_rna_mat[,j])*pred_diff_mat[,j]
}
colnames(atac_mat3) <- colnames(rna_mat2)

which(colnames(rna_mat2) == "EOMES")
round(cbind(pred_rna_mat[,"EOMES"], atac_mat3[,"EOMES"])[round(seq(1,n,length.out=100)),],2)

#################

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
color_vec <- rev(grDevices::colorRampPalette(base_palette)(n))

for(k in 1:3){
  genes <- mentioned_genes[((k-1)*20+1):min((k*20), length(mentioned_genes))]
  
  png(paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_greenleaf_development_leafplots_v3_", 
             k, ".png"),
      height = 3500, width = 3000, units = "px", res = 300)
  par(mfrow = c(5,4), mar = c(4,4,4,0.5))
  for(gene in genes){
    set.seed(10)
    idx <- sample(1:n)
    
    x_vec1 <- atac_mat3[,gene]
    y_vec1 <- pred_rna_mat[,gene]
    
    x_vec1 <- scale(x_vec1)
    y_vec1 <- scale(y_vec1)
    
    plot(x = x_vec1[idx], y = y_vec1[idx],
         pch = 16, col = color_vec[idx], main = gene, cex = 0.5,
         xlim = quantile(x_vec1, probs = c(0.01, 0.99)),
         ylim = quantile(y_vec1, probs = c(0.01, 0.99)),
         xlab = "ATAC", ylab = "RNA",
         typ = "n")
    
    points(x = x_vec1[idx], 
           y = y_vec1[idx], 
           pch = 16, col = "white", cex = 2)
    points(x = x_vec1[idx], 
           y = y_vec1[idx], 
           pch = 16, col = color_vec[idx], cex = 1.5)
    
    points(x = x_vec1[1], 
           y = y_vec1[1],
           pch = 16, col = "white", cex = 3)
    points(x = x_vec1[1], 
           y = y_vec1[1], 
           pch = 16, col = color_vec[1], cex = 2.5)
    
    x_right <- mean(x_vec1[x_vec1 >= quantile(x_vec1, probs = 0.9)])
    x_left <- mean(x_vec1[x_vec1 <= quantile(x_vec1, probs = 0.1)])
    y_top <- mean(y_vec1[y_vec1 >= quantile(y_vec1, probs = 0.9)])
    y_bot <- mean(y_vec1[y_vec1 <= quantile(y_vec1, probs = 0.1)])
    lines(c(2*x_right,2*x_left), c(2*y_top, 2*y_bot), lwd = 2, lty = 2)
  }
  graphics.off()
}
