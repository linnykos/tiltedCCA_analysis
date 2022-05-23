rm(list=ls())
load("../../../out/main/10x_mouseembryo_tiltedcca.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

rna_common <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(rna = rna_common[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec^val, rank(alignment_vec))
})
mbrain$alignment <- alignment_vec^(scaling_grid[which.max(scaling_quality)])
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(mbrain$alignment), max(mbrain$alignment), length.out = num_color)
color_vec <- sapply(mbrain$alignment, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})

png(paste0("../../../out/figures/main/10x_mouseembryo_tcca_steadystate-full.png"),
    height = 3000, width = 3000, units = "px", res = 500)
par(mar = c(4,4,4,0.5))
plot(x = mbrain[["common_tcca"]]@cell.embeddings[,1],
     y = mbrain[["common_tcca"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     xlab = colnames(mbrain[["common_tcca"]]@cell.embeddings)[1],
     ylab = colnames(mbrain[["common_tcca"]]@cell.embeddings)[2],
     main = paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nAlignment between ATAC and common RNA"),
     xaxt = "n", yaxt = "n", bty = "n")
axis(side = 1)
axis(side = 2)
graphics.off()

###########

Seurat::DefaultAssay(mbrain) <- "SCT"
mat_1 <- Matrix::t(mbrain[["SCT"]]@data[Seurat::VariableFeatures(object = mbrain),])
mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}
alignment_vec_alt <- sapply(1:n, function(i){
  df <- data.frame(rna = mat_1b[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec_alt^val, rank(alignment_vec_alt))
})
mbrain$alignment_alt <- alignment_vec_alt^(scaling_grid[which.max(scaling_quality)])
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(mbrain$alignment_alt), max(mbrain$alignment_alt), length.out = num_color)
color_vec <- sapply(mbrain$alignment_alt, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})

png(paste0("../../../out/figures/main/10x_mouseembryo_tcca_steadystate_alt-full.png"),
    height = 3000, width = 3000, units = "px", res = 500)
par(mar = c(4,4,4,0.5))
plot(x = mbrain[["common_tcca"]]@cell.embeddings[,1],
     y = mbrain[["common_tcca"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     xlab = colnames(mbrain[["common_tcca"]]@cell.embeddings)[1],
     ylab = colnames(mbrain[["common_tcca"]]@cell.embeddings)[2],
     main = paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nAlignment between ATAC and RNA"),
     xaxt = "n", yaxt = "n", bty = "n")
axis(side = 1)
axis(side = 2)
graphics.off()