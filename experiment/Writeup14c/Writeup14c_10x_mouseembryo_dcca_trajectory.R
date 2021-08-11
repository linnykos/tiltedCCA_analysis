rm(list=ls()); 

library(Seurat); library(Signac); library(multiomicCCA)
load("../../../../out/Writeup14c/10x_mouseembryo_dcca.RData")

set.seed(10)

gene_val <- sapply(gene_smoothed, function(x){
  (x$d_variance - x$c_variance)/x$e_variance
})
threshold <- 0.6
gene_bool <- sapply(gene_smoothed, function(x){min(x$c_r2, x$d_r2) > threshold})

gene_common_idx <- intersect(which(gene_bool), which(gene_val < -0.15))
gene_distinct_idx <- intersect(which(gene_bool), which(gene_val > 0.15))
length(gene_common_idx)
length(gene_distinct_idx)

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, verbose = T)
mat_1_denoised <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
rm(list = c("dcca_decomp"))

c_eig <- mbrain2[["clap"]]@cell.embeddings
d_eig <- mbrain2[["dlap"]]@cell.embeddings

tmp <- sapply(gene_common_idx, function(j){
  c_res <- compute_smooth_signal(mat_1_denoised[,j], c_eig)
  min_val <- min(c_res$smoothed_vec)
  max_val <- max(c_res$smoothed_vec)
  
  (c_res$smoothed_vec-min_val)/(max_val - min_val)
})
gene_common_vec <- rowSums(tmp)

tmp <- sapply(gene_distinct_idx, function(j){
  d_res <- compute_smooth_signal(mat_1_denoised[,j], d_eig)
  min_val <- min(d_res$smoothed_vec)
  max_val <- max(d_res$smoothed_vec)
  (d_res$smoothed_vec-min_val)/(max_val - min_val)
})
gene_distinct_vec <- rowSums(tmp)

tmp_mat <- cbind(gene_common_vec, gene_distinct_vec)
rownames(tmp_mat) <- rownames(mbrain2@meta.data)
colnames(tmp_mat) <- paste0("gene_", 1:2)

mbrain2[["gene"]] <- Seurat::CreateDimReducObject(embedding = tmp_mat, key = "gene_", assay = "RNA")

plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("gene_", 1:2), reduction = "both",
                             combine = F)
plot1[[1]] <- plot1[[1]] + ggplot2::ggtitle("Total high-common genes")
plot1[[2]] <- plot1[[2]] + ggplot2::ggtitle("Total high-distinct genes")
for(i in 1:2){
  plot1[[i]] <-plot1[[i]] + ggplot2::theme_classic()
}
tmp2 <- cowplot::plot_grid(plot1[[1]], plot1[[2]])
cowplot::save_plot(filename = "../../../../out/figures/Writeup14c/Writeup14c_10x_mouseembryo_rna_trajectory_umap.png",
                   tmp2, ncol = 2, nrow = 1, base_height = 3.5, base_asp = 4/3, device = "png")



