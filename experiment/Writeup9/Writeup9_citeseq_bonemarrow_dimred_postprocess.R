rm(list=ls())
load("../../out/Writeup9_citeseq_bonemarrow_dimred_onlyembedding.RData")

dim(umap_clr_adt)
dim(umap_log_rna)
dim(pca_clr_adt)
dim(pca_log_rna)

set.seed(10)
clust_res <- stats::kmeans(umap_log_rna, centers = 12)

par(mfrow = c(1,2))
plot(umap_log_rna[,1], umap_log_rna[,2], asp = T, pch = 16, col = clust_res$cluster)
plot(umap_clr_adt[,1], umap_clr_adt[,2], asp = T, pch = 16, col = clust_res$cluster)

par(mfrow = c(1,2))
plot(umap_log_rna[,1], umap_log_rna[,2], asp = T, pch = 16, col = clust_res$cluster)
plot(pca_log_rna[,1], pca_log_rna[,2], asp = T, pch = 16, col = clust_res$cluster)

par(mfrow = c(1,2))
plot(pca_log_rna[,1], pca_log_rna[,3], asp = T, pch = 16, col = clust_res$cluster)
plot(pca_log_rna[,1], pca_log_rna[,4], asp = T, pch = 16, col = clust_res$cluster)

par(mfrow = c(1,2))
plot(umap_clr_adt[,1], umap_clr_adt[,2], asp = T, pch = 16, col = clust_res$cluster)
plot(pca_clr_adt[,1], pca_clr_adt[,2], asp = T, pch = 16, col = clust_res$cluster)

par(mfrow = c(1,2))
plot(pca_clr_adt[,1], pca_clr_adt[,3], asp = T, pch = 16, col = clust_res$cluster)
plot(pca_clr_adt[,1], pca_clr_adt[,4], asp = T, pch = 16, col = clust_res$cluster)
