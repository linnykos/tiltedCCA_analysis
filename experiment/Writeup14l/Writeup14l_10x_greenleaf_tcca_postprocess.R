rm(list=ls())
load("../../../../out/Writeup14l/Writeup14l_10x_greenleaf_tcca.RData")

library(Seurat)

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                        group.by = "celltype", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf_umap_common-tcca.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <-Seurat::DimPlot(greenleaf, reduction = "distinct1_tcca",
                        group.by = "celltype", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nDistinct 1 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf_umap_distinct1-tcca.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <-Seurat::DimPlot(greenleaf, reduction = "distinct2_tcca",
                        group.by = "celltype", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nDistinct 2 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf_umap_distinct2-tcca.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

###############

rna_common <- multiSVD_obj$common_mat_1
rna_distinct <- multiSVD_obj$distinct_mat_1
n <- nrow(rna_common)

alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(common = rna_common[i,],
                   everything = rna_common[i,]+rna_distinct[i,])
  lm_res <- stats::lm(everything ~ common - 1, data = df)
  summary(lm_res)$r.squared
})
quantile(alignment_vec)

greenleaf$alignment <- alignment_vec^6
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "umap.wnn")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCell's RNA alignment"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$alignment <- rank(alignment_vec)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "umap.wnn")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCell's RNA alignment (Rank)"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf_umap_rna_cell_alignment.png"),
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
                            reduction = "umap.wnn")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCell's SNNs Jaccard"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

greenleaf$jaccard <- rank(cell_jaccard)
plot2 <-Seurat::FeaturePlot(greenleaf, feature = "jaccard",
                            reduction = "umap.wnn")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCell's SNNs Jaccard (Rank)"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf_umap_rna_cell_snnIntersection.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")