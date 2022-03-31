rm(list=ls())
load("../../../../out/Writeup14l/Writeup14l_10x_mouseembryo_subset_tcca.RData")
library(Seurat); library(Signac)

plot1 <-Seurat::DimPlot(mbrain, reduction = "common_tcca",
                        group.by = "label_Savercat", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_umap_common-tcca.png"),
                plot1, device = "png", width = 5.5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(mbrain, reduction = "distinct1_tcca",
                        group.by = "label_Savercat", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nDistinct 1 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_umap_distinct1-tcca.png"),
                plot1, device = "png", width = 5.5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(mbrain, reduction = "distinct2_tcca",
                        group.by = "label_Savercat", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nDistinct 2 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_umap_distinct2-tcca.png"),
                plot1, device = "png", width = 5.5, height = 5, units = "in")

################################

plot1 <-Seurat::DimPlot(mbrain, reduction = "common_tcca",
                        group.by = "Phase", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nCell-cycle phase"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_umap_cellcycle.png"),
                plot1, device = "png", width = 5.5, height = 5, units = "in")


plot1 <- Seurat::DimPlot(mbrain, reduction = "common_tcca",
                         cells.highlight = colnames(mbrain)[which(mbrain$label_Savercat == "Radial glia")])
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRadial glia")) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot2 <- Seurat::DimPlot(mbrain, reduction = "common_tcca",
                         cells.highlight = colnames(mbrain)[which(mbrain$label_Savercat == "Glioblast")])
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nGlioblast")) + Seurat::NoLegend()
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot3 <- Seurat::DimPlot(mbrain, reduction = "common_tcca",
                         cells.highlight = colnames(mbrain)[which(mbrain$label_Savercat == "Oligodendrocyte")])
plot3 <- plot3 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nOligodendrocyte")) + Seurat::NoLegend()
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot4 <- Seurat::DimPlot(mbrain, reduction = "common_tcca",
                         cells.highlight = colnames(mbrain)[which(mbrain$label_Savercat == "Neuroblast")])
plot4 <- plot4 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nNeuroblast")) + Seurat::NoLegend()
plot4 <- plot4 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot5 <- Seurat::DimPlot(mbrain, reduction = "common_tcca",
                         cells.highlight = colnames(mbrain)[which(mbrain$label_Savercat == "Forebrain GABAergic")])
plot5 <- plot5 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nForebrain GABAergic")) + Seurat::NoLegend()
plot5 <- plot5 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot6 <- Seurat::DimPlot(mbrain, reduction = "common_tcca",
                         cells.highlight = colnames(mbrain)[which(mbrain$label_Savercat == "Cortical or hippocampal glutamatergic")])
plot6 <- plot6 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCortical or hippocampal glutamatergic")) + Seurat::NoLegend()
plot6 <- plot6 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot7 <- cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 3)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_umap_celltypes.png"),
                plot7, device = "png", width = 15, height = 10, units = "in")


################################

rna_common <- multiSVD_obj$common_mat_1
rna_distinct <- multiSVD_obj$distinct_mat_1
n <- nrow(rna_common)
Seurat::DefaultAssay(mbrain) <- "SCT"
mat_1 <- Matrix::t(mbrain[["SCT"]]@data[Seurat::VariableFeatures(object = mbrain),])
mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}
mat_1b <- scale(mat_1b)

alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(common = rna_common[i,],
                   everything = mat_1b[i,])
  lm_res <- stats::lm(everything ~ common, data = df)
  summary(lm_res)$r.squared
})
quantile(alignment_vec)

mbrain$alignment <- alignment_vec
plot1 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nCell's RNA alignment"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

mbrain$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nCell's RNA alignment (Rank)"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_umap_rna_cell_alignment.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

################################

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
png("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_gene_smoothness.png",
    height = 2500, width = 2500, res = 300, units = "px")
plot(gene_var, gene_alignment, 
     main = "Mouse Embryo (10x, RNA+ATAC)\nGene smoothness vs. SD",
     xlab = "Gene stadard deviation",
     ylab = "Gene's smoothness wrt common SNN",
     pch = 16, col = rgb(0.5,0.5,0.5,0.5))
graphics.off()

seurat_tmp <- mbrain
names(gene_alignment) <- colnames(mat_1)
gene_vec <- names(gene_alignment)[order(gene_alignment, decreasing = T)[1:5]]
for(gene_name in gene_vec){
  tmp1 <- mat_1[,gene_name]
  tmp2 <- rna_common[,gene_name]
  seurat_tmp$tmp1 <- tmp1
  seurat_tmp$tmp2 <- tmp2
  plot1 <-Seurat::FeaturePlot(seurat_tmp, feature = "tmp1",
                              reduction = "umap.wnn")
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nGene: ", gene_name))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  
  plot2 <-Seurat::FeaturePlot(seurat_tmp, feature = "tmp2",
                              reduction = "umap.wnn")
  plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nGene's Common: ", gene_name))
  plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  
  plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_gene_highsnn_", gene_name, ".png"),
                  plot3, device = "png", width = 10, height = 5, units = "in")
}


seurat_tmp <- mbrain
names(gene_alignment) <- colnames(mat_1)
gene_vec <- names(gene_alignment)[order(gene_alignment, decreasing = F)[1:5]]
for(gene_name in gene_vec){
  tmp1 <- mat_1[,gene_name]
  tmp2 <- rna_common[,gene_name]
  seurat_tmp$tmp1 <- tmp1
  seurat_tmp$tmp2 <- tmp2
  plot1 <-Seurat::FeaturePlot(seurat_tmp, feature = "tmp1",
                              reduction = "umap.wnn")
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nGene: ", gene_name))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  
  plot2 <-Seurat::FeaturePlot(seurat_tmp, feature = "tmp2",
                              reduction = "umap.wnn")
  plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nGene's Common: ", gene_name))
  plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  
  plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_gene_lowsnn_", gene_name, ".png"),
                  plot3, device = "png", width = 10, height = 5, units = "in")
}

gene_vec <- names(gene_var)[which(gene_var >= 0.5)]

alignment_vec2 <- sapply(1:n, function(i){
  df <- data.frame(common = rna_common[i,gene_vec],
                   everything = rna_common[i,gene_vec]+rna_distinct[i,gene_vec])
  lm_res <- stats::lm(everything ~ common - 1, data = df)
  summary(lm_res)$r.squared
})
quantile(alignment_vec2)


mbrain$alignment <- alignment_vec2^6
plot1 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "umap.wnn")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nCell's RNA alignment (gene subset)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_umap_rna_cell_alignment_genesubset.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

##################################

snn_1 <- multiSVD_obj$snn_list$snn_1
snn_2 <- multiSVD_obj$snn_list$snn_2
n <- nrow(snn_1)
cell_jaccard <- sapply(1:n, function(i){
  idx1 <- tiltedCCA:::.nonzero_col(snn_1, col_idx = i, bool_value = F)
  idx2 <- tiltedCCA:::.nonzero_col(snn_2, col_idx = i, bool_value = F)
  
  length(intersect(idx1, idx2))/length(unique(c(idx1, idx2)))
})
quantile(cell_jaccard)

mbrain$jaccard <- cell_jaccard
plot1 <-Seurat::FeaturePlot(mbrain, feature = "jaccard",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nCell's SNNs Jaccard"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

mbrain$jaccard <- rank(cell_jaccard)
plot2 <-Seurat::FeaturePlot(mbrain, feature = "jaccard",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nCell's SNNs Jaccard (Rank)"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_umap_rna_cell_snnIntersection.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

##################################

Seurat::DefaultAssay(mbrain) <- "SCT"
mat_1 <- Matrix::t(mbrain[["SCT"]]@data[Seurat::VariableFeatures(object = mbrain),])
mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}
mat_1b <- scale(mat_1b)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)

atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
n <- nrow(mat_1b)
alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(rna = mat_1b[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})
quantile(alignment_vec)

mbrain$alignment <- alignment_vec
plot1 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\n(Everything ATAC)-(everything RNA)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

mbrain$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_everythingAtacEverythingRnaAlignment.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

##################################

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

mbrain$alignment <- alignment_vec^6
plot1 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\n(Everything ATAC)-(Common RNA)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

mbrain$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_everythingAtacCommonRnaAlignment.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")


##############################

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


mbrain$alignment <- alignment_vec
plot1 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nCommon to Distinct RNAs' correlation"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

mbrain$alignment <- rank(alignment_vec_smoothed)
plot2 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_commonCorrelationDistinct.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")


mbrain$alignment <- alignment_vec_smoothed
plot1 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nCommon-Distinct RNAs' corr. (Smoothed)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

mbrain$alignment <- rank(alignment_vec_smoothed)
plot2 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_commonCorrelationDistinct_smoothed.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

###########################

rna_common <- multiSVD_obj$common_mat_1
rna_distinct <- multiSVD_obj$distinct_mat_1
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(common = rna_common[i,],
                   distinct = rna_distinct[i,])
  lm_res <- stats::lm(distinct ~ common, data = df)
  1-tiltedCCA:::.l2norm(lm_res$residuals)/tiltedCCA:::.l2norm(rna_distinct[i,])
})
quantile(alignment_vec)

alignment_vec_smoothed <- sapply(1:n, function(i){
  idx <- c(tiltedCCA:::.nonzero_col(snn_mat, col_idx = i, bool_value = F), i)
  mean(alignment_vec[idx])
})

mbrain$alignment <- alignment_vec_smoothed
plot1 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nDistinct resid. (predict using Common, Smoothed)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

mbrain$alignment <- rank(alignment_vec_smoothed)
plot2 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo (10x, RNA+ATAC)\nRank"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_mouseembryo_subset_RNADistinctResidualsByRNACommon.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")



