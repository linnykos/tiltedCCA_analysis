rm(list=ls()); set.seed(10); gcinfo(TRUE)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

date_of_run <- Sys.time(); session_info <- sessionInfo()
load("../../out/Writeup10_10x_pbmc_preprocess4.RData"); pbmc[["ATAC"]] = NULL
load("../../out/Writeup10_10x_pbmc_dcca.RData")

##########################

res <- dcca_decomposition(dcca_res, rank_12 = K)

tmp <- dcca_res$common_factors
rownames(tmp) <- colnames(pbmc)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCA_common_')
pbmc[["common_factor"]] <- zz

png("../../out/Writeup10_pbmc_10x_dcca_common_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Common UMAP (D-CCA)")
graphics.off()

###

tmp <- RSpectra::svds(res$common_mat_1 + res$distinct_mat_1, k = K)
tmp <- multiomicCCA:::.mult_mat_vec(tmp$u, tmp$d)
rownames(tmp) <- colnames(pbmc)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCA_common_')
pbmc[["common_factor"]] <- zz

png("../../out/Writeup10_pbmc_10x_dcca_rna_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA-only UMAP (D-CCA)")
graphics.off()

###

tmp <- RSpectra::svds(res$common_mat_2 + res$distinct_mat_2, k = K)
tmp <- multiomicCCA:::.mult_mat_vec(tmp$u, tmp$d)
rownames(tmp) <- colnames(pbmc)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCA_common_')
pbmc[["common_factor"]] <- zz

png("../../out/Writeup10_pbmc_10x_dcca_atac_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("ATAC-only UMAP (D-CCA)")
graphics.off()

####

tmp <- RSpectra::svds(res$distinct_mat_1, k = K)
tmp <- multiomicCCA:::.mult_mat_vec(tmp$u, tmp$d)
rownames(tmp) <- colnames(pbmc)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCA_common_')
pbmc[["common_factor"]] <- zz

png("../../out/Writeup10_pbmc_10x_dcca_rna_distinct_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA distinct UMAP (D-CCA)")
graphics.off()

####

tmp <- RSpectra::svds(res$distinct_mat_2, k = K)
tmp <- multiomicCCA:::.mult_mat_vec(tmp$u, tmp$d)
rownames(tmp) <- colnames(pbmc)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCA_common_')
pbmc[["common_factor"]] <- zz

png("../../out/Writeup10_pbmc_10x_dcca_atac_distinct_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("ATAC distinct UMAP (D-CCA)")
graphics.off()

####

tmp <- vector("list", length = 4)
tmp[[1]] <- RSpectra::svds(res$common_mat_1, k = K)
tmp[[2]] <- RSpectra::svds(res$common_mat_2, k = K)
tmp[[3]] <- RSpectra::svds(res$distinct_mat_1, k = K)
tmp[[4]] <- RSpectra::svds(res$distinct_mat_2, k = K)

tmp <- do.call(cbind, lapply(tmp, function(svd_tmp){
  multiomicCCA:::.mult_mat_vec(svd_tmp$u, svd_tmp$d)
}))
rownames(tmp) <- colnames(pbmc)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCA_common_')
pbmc[["common_factor"]] <- zz

png("../../out/Writeup10_pbmc_10x_dcca_everything_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Everything UMAP (D-CCA)")
graphics.off()

##################################



