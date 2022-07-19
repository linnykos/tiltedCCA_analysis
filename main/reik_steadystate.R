rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_reik_tiltedcca.RData")
source("reik_colorPalette.R")

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
round(quantile(alignment_vec), 2)

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec^val, rank(alignment_vec))
})
reik$alignment <- alignment_vec^(scaling_grid[which.max(scaling_quality)])
reik$alignment2 <- alignment_vec
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(reik$alignment), max(reik$alignment), length.out = num_color)
color_vec <- sapply(reik$alignment, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})

save(color_vec, alignment_vec,
     scaling_grid, scaling_quality,
     date_of_run, session_info,
     file = "../../../out/main/10x_reik_steadystate.RData")

png(paste0("../../../out/figures/main/10x_reik_tcca_steadystate-full.png"),
    height = 3000, width = 3000, units = "px", res = 500)
par(mar = c(4,4,4,0.5))
plot(x = reik[["common_tcca"]]@cell.embeddings[,1],
     y = reik[["common_tcca"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     xlab = colnames(reik[["common_tcca"]]@cell.embeddings)[1],
     ylab = colnames(reik[["common_tcca"]]@cell.embeddings)[2],
     main = paste0("Mouse embryo (10x, RNA+ATAC)\nAlignment between ATAC and common RNA"),
     xaxt = "n", yaxt = "n", bty = "n")
axis(side = 1)
axis(side = 2)
graphics.off()

plot1 <- Seurat::VlnPlot(reik, features = "alignment2", 
                         group.by = 'celltype', 
                         sort = TRUE, pt.size = 0.1, cols = col_palette) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_steadystate_violinplot.png"),
                plot1, device = "png", width = 20, height = 5, units = "in")

