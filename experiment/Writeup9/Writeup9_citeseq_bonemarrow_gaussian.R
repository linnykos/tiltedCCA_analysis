rm(list=ls())
set.seed(10)

date_of_run <- Sys.time()
session_info <- sessionInfo()

library(Seurat)
load("../../out/Writeup9_citeseq_bonemarrow_dimred_all.RData")

mat <- t(bm[["RNA"]]@scale.data)
dim(mat)
K <- 8
bm2 <- bm

###############################################

# first try the normal PCA
svd_res <- RSpectra::svds(mat, k = K)
pca_res <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
rownames(pca_res) <- colnames(bm)

set.seed(10)
zz <- Seurat::RunUMAP(pca_res, reduction.key = 'PCAumap_')
bm2[["custom_pca"]] <- zz
png("../../out/Writeup9_citeseq_bonemarrow_pca_custom_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm2, reduction = 'custom_pca', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + NoLegend()
plot1 + ggplot2::ggtitle("PCA with 8 dimensions")
graphics.off()

# try ePCA
# we need to get the unscaled but filtered data
mat2 <- t(as.matrix(bm[["RNA"]]@data[bm[["RNA"]]@var.features,]))

dim(mat2)
epca_eig <- multiomicCCA::epca(mat2, K)
epca_res <- mat2 %*% epca_eig
set.seed(10)
zz <- Seurat::RunUMAP(epca_res, reduction.key = 'ePCAumap_')
bm2[["epca"]] <- zz
png("../../out/Writeup9_citeseq_bonemarrow_epca_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm2, reduction = 'epca', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + NoLegend()
plot1 + ggplot2::ggtitle("ePCA with 8 dimensions")
graphics.off()

# try heteroPCA
heteroPCA_cov <- multiomicCCA::heteroPCA(mat2, K)
heteroPCA_eig <- RSpectra::svds(heteroPCA_cov, k = K)
heteroPCA_res <- mat2 %*% heteroPCA_eig$u
set.seed(10)
zz <- Seurat::RunUMAP(heteroPCA_res, reduction.key = 'heteroPCAumap_')
bm2[["heteroPCA"]] <- zz
png("../../out/Writeup9_citeseq_bonemarrow_heteroPCA_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm2, reduction = 'heteroPCA', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + NoLegend()
plot1 + ggplot2::ggtitle("heteroPCA with 8 dimensions")
graphics.off()

########################################

set.seed(10)
diagnostic_res <- multiomicCCA:::softImpute_diagnostic(mat, K = 400, lambda = 0.01)
sig <- stats::sd(diagnostic_res$testing_mat[,1] - diagnostic_res$testing_mat[,2])
png("../../out/Writeup9_citeseq_bonemarrow_diagnostic.png", height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,2))
multiomicCCA::plot_prediction_against_observed(diagnostic_res$training_mat, main = "Training data",
                                               transparency = 0.05, scalar = sig, xlim = c(-5,10), ylim = c(-5,10))
multiomicCCA::plot_prediction_against_observed(diagnostic_res$testing_mat, main = "Testing data",
                                               transparency = 0.05, scalar = sig, xlim = c(-5,10), ylim = c(-5,10))
graphics.off()
