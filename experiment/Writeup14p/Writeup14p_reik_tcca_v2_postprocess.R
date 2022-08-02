rm(list=ls())
library(Seurat)

load("../../../../out/Writeup14p/Writeup14p_10x_reik_tiltedcca_v2.RData")
source("../../main/reik_colorPalette.R")

plot1 <- Seurat::DimPlot(reik, reduction = "common_tcca",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_reik_tcca_v2-umap_common.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")

plot1 <- Seurat::DimPlot(reik, reduction = "distinct1_tcca",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nRNA Distinct subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_reik_tcca_v2-umap_distinct1.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")

plot1 <- Seurat::DimPlot(reik, reduction = "distinct2_tcca",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nATAC Distinct subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_reik_tcca_v2-umap_distinct2.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")

################3######################


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

png(paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_reik_tcca_steadystate-full.png"),
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
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_reik_steadystate_violinplot.png"),
                plot1, device = "png", width = 20, height = 5, units = "in")



