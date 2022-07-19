rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/main/10x_greenleaf_tcca_RNA-ATAC.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

svd_1 <- tiltedCCA:::.svd_safe(mat = multiSVD_obj$common_mat_1,
                               check_stability = T,
                               K = ncol(multiSVD_obj$svd_2$u), # positive integer
                               mean_vec = NULL, # boolean, NULL or vector
                               rescale = F, # boolean
                               scale_max = NULL, # NULL or positive integer
                               sd_vec = NULL)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
rna_everything <- multiSVD_obj$common_mat_1 + multiSVD_obj$distinct_mat_1

n <- nrow(rna_everything)
alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(rna = rna_everything[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec^val, rank(alignment_vec))
})
alignment_vec2 <- alignment_vec^(scaling_grid[which.max(scaling_quality)])
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(alignment_vec2), max(alignment_vec2), length.out = num_color)
color_vec <- sapply(alignment_vec2, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})

png(paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_greenleaf_tcca_steadystate_v2.png"),
    height = 3000, width = 3000, units = "px", res = 500)
par(mar = c(4,4,4,0.5))
plot(x = greenleaf[["common_tcca"]]@cell.embeddings[,1],
     y = greenleaf[["common_tcca"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     xlab = colnames(greenleaf[["common_tcca"]]@cell.embeddings)[1],
     ylab = colnames(greenleaf[["common_tcca"]]@cell.embeddings)[2],
     main = paste0("Human brain (10x, RNA+ATAC)\nAlignment b/w ATAC (rotated to common RNA) and RNA"),
     xaxt = "n", yaxt = "n", bty = "n")
axis(side = 1)
axis(side = 2)
graphics.off()