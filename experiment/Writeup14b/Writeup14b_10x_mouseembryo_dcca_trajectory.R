rm(list=ls()); 

library(Seurat); library(Signac); library(multiomicCCA)
load("../../../../out/Writeup14b/Writeup14b_10x_mouseembryo_dcca_laplacian_variablecalculations.RData")

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

load("../../../../out/Writeup14b/Writeup14b_10x_mouseembryo_dcca_laplacian_variablecalculations_tmp.RData")
