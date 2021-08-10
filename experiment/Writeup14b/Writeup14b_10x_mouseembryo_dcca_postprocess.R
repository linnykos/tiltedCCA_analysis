rm(list=ls())
load("../../../../out/Writeup14b/Writeup14b_10x_mouseembryo_dcca_laplacian_variablecalculations.RData")

val_vec <- sapply(gene_smoothed, function(x){
  (x$d_variance - x$c_variance)/x$e_variance
})
name_vec <- rownames(dcca_res$svd_1$v)
p1 <- length(name_vec)
factor_vec <- rep(0, p1); 
factor_vec[which(sapply(gene_smoothed, function(x){min(x$c_r2, x$d_r2) > 0.6}))] <- 2
idx <- which(factor_vec == 2)
idx2 <- which.max(val_vec[idx])
factor_vec[idx[idx2]] <- 1; factor_vec <- as.factor(factor_vec)
col_vec <- c("black", "red", "green"); names(col_vec) <- c("0", "1", "2")
tmp <- multiomicCCA::plot_laplacian_variables(val_vec, name_vec, factor_vec, col_vec,
                                              ylab = "Distinct-common (normalized)", main = "RNA enrichment",
                                              text_cex = 4)
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_rna_var_enrichment.png",
                tmp, device = "png", width = 3, height = 3, units = "in")

################

val_vec <- sapply(atac_smoothed, function(x){
  (x$d_variance - x$c_variance)/x$e_variance
})
name_vec <- rownames(dcca_res$svd_2$v)
p2 <- length(name_vec)
factor_vec <- rep(0, p2); 
factor_vec[which(sapply(atac_smoothed, function(x){min(x$c_r2, x$d_r2) > 0.6}))] <- 2
idx <- which(factor_vec == 2)
idx2 <- which.max(val_vec[idx])
factor_vec[idx[idx2]] <- 1; factor_vec <- as.factor(factor_vec)
col_vec <- c("black", "red", "green"); names(col_vec) <- c("0", "1", "2")
tmp <- multiomicCCA::plot_laplacian_variables(val_vec, name_vec, factor_vec, col_vec,
                                              ylab = "Distinct-common (normalized)", main = "ATAC enrichment",
                                              text_cex = 4)
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_atac_var_enrichment.png",
                tmp, device = "png", width = 3, height = 3, units = "in")







