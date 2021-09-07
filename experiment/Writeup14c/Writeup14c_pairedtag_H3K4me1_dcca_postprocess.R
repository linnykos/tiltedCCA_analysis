rm(list=ls())
load("../../../../out/Writeup14c/Writeup14c_pairedtag_H3K4me1_dcca.RData")

title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")

# plot RNA embeddings
for(i in 1:3){
  plot1 <- Seurat::DimPlot(pairedtag2, reduction = main_vec[i],
                           group.by = "celltype", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag H3K4me1 (RNA)\n", title_vec[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14c/Writeup14c_pairedtag_H3K4me1_dcca_rna_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# plot ATAC embeddings
for(i in 1:3){
  plot1 <- Seurat::DimPlot(pairedtag2, reduction = paste0(main_vec[i], "2"),
                           group.by = "celltype", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag H3K4me1 (Protein)\n", title_vec[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14c/Writeup14c_pairedtag_H3K4me1_dcca_histone_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# plot the combined view
plot1 <- Seurat::DimPlot(pairedtag2, reduction = "combined",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Paired-Tag H3K4me1\nBoth Common")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14c/Writeup14c_pairedtag_H3K4me1_dcca_both_common_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

# plot the everything combined for both
plot1 <- Seurat::DimPlot(pairedtag2, reduction = "both",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Paired-Tag H3K4me1\nBoth Everything")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14c/Writeup14c_pairedtag_H3K4me1_dcca_both_everything_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


#######################################
# now plot the bases -- first RNA
plot1 <- Seurat::FeaturePlot(pairedtag2, features = paste0("clap_", 1:16), reduction = "combined")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14c/Writeup14c_pairedtag_H3K4me1_dcca_rna_common_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(pairedtag2, features = paste0("dlap_", 1:16), reduction = "distinct")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14c/Writeup14c_pairedtag_H3K4me1_dcca_rna_distinct_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(pairedtag2, features = paste0("elap_", 1:16), reduction = "both")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14c/Writeup14c_pairedtag_H3K4me1_dcca_rna_everything_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

# next ATAC
plot1 <- Seurat::FeaturePlot(pairedtag2, features = paste0("clap2_", 1:16), reduction = "combined")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14c/Writeup14c_pairedtag_H3K4me1_dcca_atac_common_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(pairedtag2, features = paste0("dlap2_", 1:16), reduction = "distinct2")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14c/Writeup14c_pairedtag_H3K4me1_dcca_atac_distinct_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(pairedtag2, features = paste0("elap2_", 1:16), reduction = "both")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14c/Writeup14c_pairedtag_H3K4me1_dcca_atac_everything_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

##################


membership_vec <- as.factor(pairedtag2[["celltype"]][,1]); nn <- 30

# set.seed(10)
# rna_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
#                                          data_1 = T, data_2 = F,
#                                          radius_quantile = 0.5, symmetrize = F, 
#                                          bool_matrix = T, verbose = T)
# 
# set.seed(10)
# protein_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
#                                              data_1 = F, data_2 = T,
#                                              radius_quantile = 0.5,
#                                              bool_matrix = T, symmetrize = F, verbose = T)

# compute local enrichment
set.seed(10)
rna_local <- multiomicCCA::clisi_information(rna_frnn$c_g, rna_frnn$d_g, rna_frnn$e_g, 
                                             membership_vec = membership_vec)
# rna_local$common_clisi$membership_info
# rna_local$distinct_clisi$membership_info

set.seed(10)
protein_local <- multiomicCCA::clisi_information(protein_frnn$c_g, protein_frnn$d_g, protein_frnn$e_g, 
                                                 membership_vec = membership_vec)
# atac_local$common_clisi$membership_info
# atac_local$distinct_clisi$membership_info

tmp <- multiomicCCA::plot_clisi(rna_local, protein_local)
tmp2 <- cowplot::plot_grid(tmp[[1]], tmp[[2]])
cowplot::save_plot(filename = "../../../../out/figures/Writeup14c/Writeup14c_pairedtag_H3K4me1_enrichment.png", 
                   tmp2, ncol = 1, nrow = 2, base_height = 1.75, base_asp = 4, device = "png")


################################

val_vec <- sapply(gene_smoothed, function(x){
  (x$d_variance - x$c_variance)/x$e_variance
})
name_vec <- rownames(dcca_res$svd_1$v)
p1 <- length(name_vec)
factor_vec <- rep(0, p1)
threshold <- 0.6
factor_vec[which(sapply(gene_smoothed, function(x){min(x$c_r2, x$d_r2) > threshold}))] <- 2
idx <- which(factor_vec == 2)
idx2 <- which.max(val_vec[idx])
factor_vec[idx[idx2]] <- 1; factor_vec <- as.factor(factor_vec)
col_vec <- c("black", "red", "green"); names(col_vec) <- c("0", "1", "2")
tmp <- multiomicCCA::plot_laplacian_variables(val_vec, name_vec, factor_vec, col_vec,
                                              ylab = "Distinct-common (normalized)", 
                                              main = paste0("RNA enrichment\nThreshold: ", threshold),
                                              text_cex = 4)
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14c/Writeup14c_pairedtag_H3K4me1_rna_var_enrichment.png",
                tmp, device = "png", width = 3, height = 3, units = "in")

################

val_vec <- sapply(histone_smoothed, function(x){
  (x$d_variance - x$c_variance)/x$e_variance
})
name_vec <- rownames(dcca_res$svd_2$v)
p2 <- length(name_vec)
factor_vec <- rep(0, p2)
threshold <- 0.15
factor_vec[which(sapply(histone_smoothed, function(x){min(x$c_r2, x$d_r2) > threshold}))] <- 2
idx <- which(factor_vec == 2)
idx2 <- which.max(val_vec[idx])
factor_vec[idx[idx2]] <- 1; factor_vec <- as.factor(factor_vec)
col_vec <- c("black", "red", "green"); names(col_vec) <- c("0", "1", "2")
tmp <- multiomicCCA::plot_laplacian_variables(val_vec, name_vec, factor_vec, col_vec,
                                              ylab = "Distinct-common (normalized)", 
                                              main = paste0("Histone enrichment\nThreshold: ", threshold),
                                              text_cex = 4)
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14c/Writeup14c_pairedtag_H3K4me1_histone_var_enrichment.png",
                tmp, device = "png", width = 3, height = 3, units = "in")


