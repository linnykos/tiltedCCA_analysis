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

# version 1: let's do the relative difference between these two, then smoothed by npreg
diff_mat <- (rna_mat - atac_mat)/rna_mat
diff_mat <- pmin(pmax(diff_mat, -5), 5)
diff_mat <- abs(diff_mat[,mentioned_genes])

p <- ncol(diff_mat); n <- nrow(diff_mat)
pred_diff_mat <- sapply(1:p, function(j){
  print(j)
  
  tmp_df <- data.frame(y = diff_mat[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
colnames(pred_diff_mat) <- colnames(diff_mat)

for(k in 1:3){
  genes <- mentioned_genes[((k-1)*20+1):min((k*20), length(mentioned_genes))]
  
  png(paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_greenleaf_steadystate_mentionedgenes_difference_v1-",
             k, ".png"),
      height = 3500, width = 3000, units = "px", res = 300)
  par(mfrow = c(5,4), mar = c(4,4,4,0.5))
  for(gene in genes){
    plot(diff_mat[,gene], pch = 16, col = rgb(0.5, 0.5, 0.5, 0.2),
         xlab = "Pseudotime", ylab = "Relative difference",
         main = gene)
    lines(pred_diff_mat[,gene], col = "white", lwd = 3)
    lines(pred_diff_mat[,gene], col = "black", lwd = 2)
  }
  graphics.off()
}

##################################

# version 2: rescaling

rna_mat2 <- t(scale(t(rna_mat)))
atac_mat2 <- t(scale(t(atac_mat)))

diff_mat <- (rna_mat2 - atac_mat2)/rna_mat2
diff_mat <- abs(diff_mat[,mentioned_genes])
diff_mat <- pmin(diff_mat, 5)

p <- ncol(diff_mat); n <- nrow(diff_mat)
pred_diff_mat <- sapply(1:p, function(j){
  print(j)
  
  tmp_df <- data.frame(y = diff_mat[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
colnames(pred_diff_mat) <- colnames(diff_mat)

for(k in 1:3){
  genes <- mentioned_genes[((k-1)*20+1):min((k*20), length(mentioned_genes))]
  
  png(paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_greenleaf_steadystate_mentionedgenes_difference_v2-",
             k, ".png"),
      height = 3500, width = 3000, units = "px", res = 300)
  par(mfrow = c(5,4), mar = c(4,4,4,0.5))
  for(gene in genes){
    plot(diff_mat[,gene], pch = 16, col = rgb(0.5, 0.5, 0.5, 0.2),
         xlab = "Pseudotime", ylab = "Relative difference",
         main = gene)
    lines(pred_diff_mat[,gene], col = "white", lwd = 3)
    lines(pred_diff_mat[,gene], col = "black", lwd = 2)
  }
  graphics.off()
}
