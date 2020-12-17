rm(list=ls())
load("../../out/Writeup9_10x_pbmc_dimred_onlyembedding.RData")

dim(pca_log_rna)
dim(pca_sct_rna)
dim(pca_tfidf_atac)
dim(umap_sct_rna)
dim(umap_tfidf_atac)

set.seed(10)
clust_res <- stats::kmeans(umap_sct_rna, centers = 12)

par(mfrow = c(1,2))
plot(umap_sct_rna[,1], umap_sct_rna[,2], asp = T, pch = 16, col = clust_res$cluster)
plot(umap_tfidf_atac[,1], umap_tfidf_atac[,2], asp = T, pch = 16, col = clust_res$cluster)

par(mfrow = c(1,2))
plot(umap_sct_rna[,1], umap_sct_rna[,2], asp = T, pch = 16, col = clust_res$cluster)
plot(pca_sct_rna[,1], pca_sct_rna[,2], asp = T, pch = 16, col = clust_res$cluster)

par(mfrow = c(1,2))
plot(pca_sct_rna[,1], pca_sct_rna[,2], asp = T, pch = 16, col = clust_res$cluster)
plot(pca_log_rna[,1], pca_log_rna[,2], asp = T, pch = 16, col = clust_res$cluster)

par(mfrow = c(1,2))
plot(pca_log_rna[,1], pca_log_rna[,3], asp = T, pch = 16, col = clust_res$cluster)
plot(pca_log_rna[,1], pca_log_rna[,4], asp = T, pch = 16, col = clust_res$cluster)

par(mfrow = c(1,2))
plot(umap_tfidf_atac[,1], umap_tfidf_atac[,2], asp = T, pch = 16, col = clust_res$cluster)
plot(pca_tfidf_atac[,1], pca_tfidf_atac[,2], asp = T, pch = 16, col = clust_res$cluster)

par(mfrow = c(1,2))
plot(pca_tfidf_atac[,1], pca_tfidf_atac[,3], asp = T, pch = 16, col = clust_res$cluster)
plot(pca_tfidf_atac[,1], pca_tfidf_atac[,4], asp = T, pch = 16, col = clust_res$cluster)

