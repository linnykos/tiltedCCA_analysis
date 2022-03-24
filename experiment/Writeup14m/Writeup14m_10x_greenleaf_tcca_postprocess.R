rm(list=ls())
load("../../../../out/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_tcca.RData")

library(Seurat)

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_umap_common-tcca.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(greenleaf, reduction = "distinct1_tcca",
                        group.by = "celltype", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nDistinct 1 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_umap_distinct1-tcca.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(greenleaf, reduction = "distinct2_tcca",
                        group.by = "celltype", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nDistinct 2 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_umap_distinct2-tcca.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

###############

snn_1 <- multiSVD_obj$snn_list$snn_1
snn_2 <- multiSVD_obj$snn_list$snn_2
n <- nrow(snn_1)
cell_jaccard <- sapply(1:n, function(i){
  idx1 <- tiltedCCA:::.nonzero_col(snn_1, col_idx = i, bool_value = F)
  idx2 <- tiltedCCA:::.nonzero_col(snn_2, col_idx = i, bool_value = F)
  
  length(intersect(idx1, idx2))/length(unique(c(idx1, idx2)))
})

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
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_snnIntersection.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

##################################

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
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_everythingAtacEverythingRnaAlignment.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

################

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

greenleaf$alignment <- alignment_vec
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
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_everythingAtacCommonRnaAlignment.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

##########################

multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
dimred_1 <- tiltedCCA:::.get_postDimred(multiSVD_obj, averaging_mat = NULL)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
dimred_2 <- tiltedCCA:::.get_postDimred(multiSVD_obj, averaging_mat = NULL)
dimred <- cbind(dimred_1, dimred_2)
param <- multiSVD_obj$param

snn_mat <- tiltedCCA:::.form_snn_mat(mat = dimred, 
                                     num_neigh = param$snn_num_neigh,
                                     bool_cosine = param$snn_bool_cosine,
                                     bool_intersect = param$snn_bool_intersect,
                                     min_deg = param$snn_min_deg)

rna_common <- multiSVD_obj$common_mat_1
rna_distinct <- multiSVD_obj$distinct_mat_1
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  vec1 <- rna_common[i,]
  vec2 <- rna_distinct[i,]
  vec1 <- scale(vec1); vec2 <- scale(vec2)
  abs(stats::cor(vec1, vec2))
})

alignment_vec_smoothed <- sapply(1:n, function(i){
  idx <- c(tiltedCCA:::.nonzero_col(snn_mat, col_idx = i, bool_value = F), i)
  mean(alignment_vec[idx])
})

greenleaf$alignment <- alignment_vec
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCommon to Distinct RNAs' correlation"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$alignment <- rank(alignment_vec_smoothed)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_commonCorrelationDistinct.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")


greenleaf$alignment <- alignment_vec_smoothed
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCommon-Distinct RNAs' corr. (Smoothed)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$alignment <- rank(alignment_vec_smoothed)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_commonCorrelationDistinct_smoothed.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")




