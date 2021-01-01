# try 1) separate PCAs, 2) concatenate both matrices,
# 3) concatenate the principal loadings, 4) the vanilla WNN,
# 5) the vanilla D-CCA

rm(list=ls()); set.seed(10); gcinfo(TRUE)

library(Seurat); library(Signac)
load("../../out/Writeup10_10x_pbmc_preprocess4.RData")

###############

# first try separate PCAs, already calculated
png("../../out/Writeup10_pbmc_10x_seurat_rna_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "predicted.id", 
                 label = TRUE, label.size = 2.5, repel = TRUE) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA UMAP (Vanilla)")
graphics.off()

png("../../out/Writeup10_pbmc_10x_seurat_atac_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "predicted.id", 
                 label = TRUE, label.size = 2.5, repel = TRUE) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("ATAC UMAP (Vanilla)")
graphics.off()

png("../../out/Writeup10_pbmc_10x_seurat_wnn_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "predicted.id", 
                 label = TRUE, label.size = 2.5, repel = TRUE) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("WNN UMAP (Seurat)")
graphics.off()

dim(pbmc[["lsi"]])

png("../../out/Writeup10_pbmc_10x_seurat_weight.png", height = 1500, width = 2500, units = "px", res = 300)
Seurat::VlnPlot(pbmc, features = "SCT.weight", group.by = 'predicted.id', sort = TRUE, pt.size = 0.1) +
  NoLegend() + ggplot2::ggtitle("Seurat WNN's RNA weight")
graphics.off()

##################################

# pick a rank
mat_1 <- t(pbmc[["SCT"]]@scale.data)
mat_2 <- t(pbmc[["ATAC"]]@scale.data)
dim(mat_1); dim(mat_2)
svd_res_1 <- RSpectra::svds(mat_1, k = 25)
svd_res_2 <- RSpectra::svds(mat_2, k = 25)
-diff(svd_res_1$d)/svd_res_1$d[-1]; -diff(svd_res_2$d)/svd_res_2$d[-1]
K <- 7
png("../../out/Writeup10_pbmc_10x_rank.png", height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,2))
graphics::plot(svd_res_1$d, xlab = "Rank", ylab = "Singular value", main = "SVD for RNA", pch = 16)
graphics::lines(rep(K+.5, 2), c(-1e5, 1e5), col = "red", lty = 2, lwd = 2)
graphics::plot(svd_res_2$d, xlab = "Rank", ylab = "Singular value", main = "SVD for ATAC", pch = 16)
graphics::lines(rep(K+.5, 2), c(-1e5, 1e5), col = "red", lty = 2, lwd = 2)
graphics.off()

##############################
# try method 2, just concenating the matrices naively
rm(list = c("plot1", "svd_res_1", "svd_res_2")); gc(verbose = T)
svd_res_1 <- RSpectra::svds(mat_1, k = 1)
svd_res_2 <- RSpectra::svds(mat_2, k = 1)
mat_total <- cbind(mat_1/svd_res_1$d[1], mat_2/svd_res_2$d[1]) # oops! error, "cannot allocate vector of size 8.3 Gb"
pca_total <- RSpectra::svds(mat_total, k = 2*K)
tmp <- multiomicCCA:::.mult_mat_vec(pca_total$u, pca_total$d)
rownames(tmp) <- rownames(pbmc[["pca"]])
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'PCATotal_')
pbmc[["pca_total_umap"]] <- zz
png("../../out/Writeup10_pbmc_10x_pca_total_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'pca_total_umap', group.by = "predicted.id", 
                         label = TRUE, label.size = 2.5, repel = TRUE) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("PCA total UMAP")
graphics.off()

# try method 3, concatenating the low-dimen embedding
mat_1 <- pbmc[["pca"]]@cell.embeddings
mat_2 <- pbmc[["lsi"]]@cell.embeddings
svd_res_1 <- RSpectra::svds(mat_1, k = 1)
svd_res_2 <- RSpectra::svds(mat_2, k = 1)

mat_total <- cbind(mat_1/svd_res_1$d[1], mat_2/svd_res_2$d[1])
rownames(mat_total) <- rownames(pbmc[["pca"]])
set.seed(10)
zz <- Seurat::RunUMAP(mat_total, reduction.key = 'PCAConcat_')
pbmc[["pca_concate_umap"]] <- zz
png("../../out/Writeup10_pbmc_10x_pca_concate_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'pca_concate_umap', group.by = "predicted.id", 
                         label = TRUE, label.size = 2.5, repel = TRUE) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("PCA-concatenated UMAP")
graphics.off()
