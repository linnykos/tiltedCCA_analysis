rm(list=ls())
load("../../../../out/Writeup14l/Writeup14l_10x_greenleaf2_tcca.RData")

library(Seurat)

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_common-tcca.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(greenleaf, reduction = "distinct1_tcca",
                        group.by = "celltype", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nDistinct 1 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_distinct1-tcca.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(greenleaf, reduction = "distinct2_tcca",
                        group.by = "celltype", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nDistinct 2 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_distinct2-tcca.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

###############

rna_common <- multiSVD_obj$common_mat_1
rna_distinct <- multiSVD_obj$distinct_mat_1
n <- nrow(rna_common)
Seurat::DefaultAssay(greenleaf) <- "SCT"
mat_1 <- Matrix::t(greenleaf[["SCT"]]@data[Seurat::VariableFeatures(object = greenleaf),])
mat_1 <- scale(mat_1)

alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(common = rna_common[i,],
                   everything = mat_1[i,])
  lm_res <- stats::lm(everything ~ common, data = df)
  summary(lm_res)$r.squared
  
  # stats::cor(rna_common[i,], rna_common[i,]+rna_distinct[i,])
  # stats::cor(rna_common[i,], mat_1[i,])
})
quantile(alignment_vec)

greenleaf$alignment <- alignment_vec
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCell's RNA alignment"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCell's RNA alignment (Rank)"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_rna_cell_alignment.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

######################

atac_common <- multiSVD_obj$common_dimred_2
svd_2 <- multiSVD_obj$svd_2
tmp <- tiltedCCA:::.mult_mat_vec(svd_2$u, svd_2$d)
tmp <- scale(tmp)
tmp <- t(sapply(1:nrow(tmp), function(i){
  tmp[i,]/tiltedCCA:::.l2norm(tmp[i,])
}))
atac_mat <- tmp
n <- nrow(atac_common)

alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(common = atac_common[i,],
                   everything = atac_mat[i,])
  lm_res <- stats::lm(everything ~ common, data = df)
  summary(lm_res)$r.squared
})
quantile(alignment_vec)

greenleaf$alignment <- alignment_vec
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCell's ATAC alignment"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCell's ATAC alignment (Rank)"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_atac_cell_alignment.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

######################

snn_1 <- multiSVD_obj$snn_list$snn_1
snn_2 <- multiSVD_obj$snn_list$snn_2
n <- nrow(snn_1)
cell_jaccard <- sapply(1:n, function(i){
  idx1 <- tiltedCCA:::.nonzero_col(snn_1, col_idx = i, bool_value = F)
  idx2 <- tiltedCCA:::.nonzero_col(snn_2, col_idx = i, bool_value = F)
  
  length(intersect(idx1, idx2))/length(unique(c(idx1, idx2)))
})
quantile(cell_jaccard)

greenleaf$jaccard <- cell_jaccard
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "jaccard",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCell's SNNs Jaccard"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$jaccard <- rank(cell_jaccard)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "jaccard",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCell's SNNs Jaccard (Rank)"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_rna_cell_snnIntersection.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

#######################

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         group.by = "Phase", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCell-cycle phase"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_cellcycle.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

#######################

cell_idx <- which(greenleaf$celltype == "GluN4")
snn_mat <- multiSVD_obj$snn_list$snn_1
rna_nn <- sapply(1:n, function(i){
  idx <- tiltedCCA:::.nonzero_col(snn_mat, col_idx = i, bool_value = F)
  length(intersect(idx, cell_idx))
})
snn_mat <- multiSVD_obj$snn_list$snn_2
atac_nn <- sapply(1:n, function(i){
  idx <- tiltedCCA:::.nonzero_col(snn_mat, col_idx = i, bool_value = F)
  length(intersect(idx, cell_idx))
})

greenleaf$rna_nn <- rna_nn
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "rna_nn",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\n# GluN4 NN's (RNA)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$atac_nn <- rna_nn
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "atac_nn",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\n# GluN4 NN's (ATAC)"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_glun4NNs.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

## 

cell_idx <- which(greenleaf$celltype == "GluN5")
snn_mat <- multiSVD_obj$snn_list$snn_1
rna_nn <- sapply(1:n, function(i){
  idx <- tiltedCCA:::.nonzero_col(snn_mat, col_idx = i, bool_value = F)
  length(intersect(idx, cell_idx))
})
snn_mat <- multiSVD_obj$snn_list$snn_2
atac_nn <- sapply(1:n, function(i){
  idx <- tiltedCCA:::.nonzero_col(snn_mat, col_idx = i, bool_value = F)
  length(intersect(idx, cell_idx))
})

greenleaf$rna_nn <- rna_nn
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "rna_nn",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\n# GluN5 NN's (RNA)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$atac_nn <- rna_nn
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "atac_nn",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\n# GluN5 NN's (ATAC)"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_glun5NNs.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

###########################

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         cells.highlight = colnames(greenleaf)[which(greenleaf$celltype == "Cyc. Prog.")])
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCycling Progenitors")) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot2 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         cells.highlight = colnames(greenleaf)[which(greenleaf$celltype == "nIPC/GluN1")])
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nnIPC or GluN1")) + Seurat::NoLegend()
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot3 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         cells.highlight = colnames(greenleaf)[which(greenleaf$celltype == "GluN2")])
plot3 <- plot3 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nGluN2")) + Seurat::NoLegend()
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot4 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         cells.highlight = colnames(greenleaf)[which(greenleaf$celltype == "GluN3")])
plot4 <- plot4 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nGluN3")) + Seurat::NoLegend()
plot4 <- plot4 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot5 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         cells.highlight = colnames(greenleaf)[which(greenleaf$celltype == "GluN4")])
plot5 <- plot5 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nGluN4")) + Seurat::NoLegend()
plot5 <- plot5 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot6 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         cells.highlight = colnames(greenleaf)[which(greenleaf$celltype == "GluN5")])
plot6 <- plot6 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nGluN5")) + Seurat::NoLegend()
plot6 <- plot6 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot7 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         cells.highlight = colnames(greenleaf)[which(greenleaf$celltype == "RG")])
plot7 <- plot7 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRG")) + Seurat::NoLegend()
plot7 <- plot7 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot8 <- cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, ncol = 3)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_celltypes.png"),
                plot8, device = "png", width = 15, height = 15, units = "in")

#################################

common_score <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, 
                                         what = "common_score")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
score_1 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, 
                                    what = "score")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
score_2 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, 
                                    what = "score")
dimred_1 <- tcrossprod(common_score, score_1[,1:ncol(common_score)])
dimred_2 <- tcrossprod(common_score, score_2[,1:ncol(common_score)])

alignment_vec <- sapply(1:n, function(i){
  stats::cor(dimred_1[i,], dimred_2[i,])
})
quantile(alignment_vec)

greenleaf$alignment <- alignment_vec^6
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCell's ATAC-RNA common alignment"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_rnaAtacCommonAlignment.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

##########################################

# let's try something else
mat_1 <- Matrix::t(greenleaf[["SCT"]]@data[Seurat::VariableFeatures(object = greenleaf),])
mat_1 <- scale(mat_1)
common_score <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, 
                                         what = "common_score")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
score_1 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F,
                                    what = "score")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
score_2 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, 
                                    what = "score")
k <- ncol(common_score)
n <- nrow(common_score)
# rna_dimred <- common_score %*% crossprod(score_1[,1:k], tiltedCCA:::.mult_mat_vec(svd_1$u, svd_1$d))/n
# rna_dimred <- tiltedCCA:::.mult_mat_vec(svd_1$u, svd_1$d)/n
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
# atac_dimred <- common_score %*% crossprod(score_2[,1:k], tiltedCCA:::.mult_mat_vec(svd_2$u, svd_2$d)) %*% rotation_mat/n
tmp <- common_score %*% crossprod(score_2[,1:k], svd_2$u) %*% rotation_mat/n
atac_pred <-  tcrossprod(tiltedCCA:::.mult_mat_vec(tmp, svd_1$d), svd_1$v)

alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(rna = mat_1[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
  #stats::cor(rna_dimred[i,], atac_dimred[i,])
})
quantile(alignment_vec)

greenleaf$alignment <- alignment_vec
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCommon ATAC to everything RNA"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_commonAtac2everythingRna.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

############################

rna_common <- multiSVD_obj$common_mat_1
rna_distinct <- multiSVD_obj$distinct_mat_1
mat_1 <- rna_common + rna_distinct

p <- ncol(rna_common)
gene_alignment <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  df <- cbind(mat_1[,j], multiSVD_obj$laplacian_list$common_laplacian)
  df <- as.data.frame(df)
  colnames(df)[1] <- "y"
  lm_res <- stats::lm(y ~ ., data = df)
  summary(lm_res)$r.squared
})
quantile(gene_alignment)

gene_var <- matrixStats::colSds(mat_1)
names(gene_var) <- colnames(mat_1)
png("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_gene_smoothness.png",
    height = 2500, width = 2500, res = 300, units = "px")
plot(gene_var, gene_alignment, 
     main = "Human brain (10x, RNA+ATAC)\nGene smoothness vs. SD",
     xlab = "Gene stadard deviation",
     ylab = "Gene's smoothness wrt common SNN",
     pch = 16, col = rgb(0.5,0.5,0.5,0.5))
graphics.off()

######################################

Seurat::DefaultAssay(greenleaf) <- "SCT"
mat_1 <- Matrix::t(greenleaf[["SCT"]]@data[Seurat::VariableFeatures(object = greenleaf),])
mat_1 <- scale(mat_1)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)

atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
n <- nrow(mat_1)
alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(rna = mat_1[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})
quantile(alignment_vec)

greenleaf$alignment <- alignment_vec
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\n(Everything ATAC)-(everything RNA)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_everythingAtacEverythingRnaAlignment.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

#################################

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
quantile(alignment_vec)

greenleaf$alignment <- alignment_vec^6
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\n(Everything ATAC)-(Common RNA)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_everythingAtacCommonRnaAlignment.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

#########################

rna_common <- multiSVD_obj$common_mat_1
Seurat::DefaultAssay(greenleaf) <- "SCT"
mat_1 <- Matrix::t(greenleaf[["SCT"]]@data[Seurat::VariableFeatures(object = greenleaf),])
mat_1 <- scale(mat_1)
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  mean(abs(rna_common[i,]/mat_1[i,]))
})

greenleaf$alignment <- alignment_vec
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCommon to observed RNA ratio"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_commonratio.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

#####

# rna_common <- multiSVD_obj$tcca_obj$common_score
# rna_distinct <- multiSVD_obj$tcca_obj$distinct_score_1
rna_common <- multiSVD_obj$common_mat_1
rna_distinct <- multiSVD_obj$distinct_mat_1
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  vec1 <- rna_common[i,]
  vec2 <- rna_common[i,]+rna_distinct[i,]
  vec1 <- scale(vec1); vec2 <- scale(vec2)
  stats::cor(vec1, vec2)
})

greenleaf$alignment <- alignment_vec
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCommon to Common+Distinct correlation"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_commonCorrelationEverything.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

#####

# rna_common <- multiSVD_obj$tcca_obj$common_score
# rna_distinct <- multiSVD_obj$tcca_obj$distinct_score_1
rna_common <- multiSVD_obj$common_mat_1
rna_distinct <- multiSVD_obj$distinct_mat_1
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  vec1 <- rna_common[i,]
  vec2 <- rna_distinct[i,]
  vec1 <- scale(vec1); vec2 <- scale(vec2)
  abs(stats::cor(vec1, vec2))
})

greenleaf$alignment <- alignment_vec
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCommon to Distinct RNAs' correlation"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_commonCorrelationDistinct.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")


# rna_common <- multiSVD_obj$tcca_obj$common_score
# rna_distinct <- multiSVD_obj$tcca_obj$distinct_score_1
rna_common <- multiSVD_obj$tcca_obj$common_score
rna_distinct <-multiSVD_obj$tcca_obj$distinct_score_1
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  vec1 <- rna_common[i,]
  vec2 <- rna_distinct[i,]
  if(length(vec1) < length(vec2)) vec1 <- c(vec1, rep(0, length(vec2)-length(vec1)))
  vec1 <- scale(vec1); vec2 <- scale(vec2)
  abs(stats::cor(vec1, vec2))
})

greenleaf$alignment <- alignment_vec
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCommon to Distinct scores' correlation"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_umap_commonCorrelationDistinct2.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")





