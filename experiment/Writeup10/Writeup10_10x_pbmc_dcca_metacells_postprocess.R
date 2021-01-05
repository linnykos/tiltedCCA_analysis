rm(list=ls()); set.seed(10); gcinfo(TRUE)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

date_of_run <- Sys.time(); session_info <- sessionInfo()
load("../../out/Writeup10_10x_pbmc_preprocess4.RData"); pbmc[["ATAC"]] = NULL
load("../../out/Writeup10_10x_pbmc_dcca_metacells2.RData")

#############

res <- dcca_decomposition(dcca_res, rank_c = K)
rownames(res$common_mat_1) <- colnames(pbmc)

set.seed(10)
tmp <- extract_embedding(res, common = T)
tmp <- dcca_res$common_factors
rownames(tmp) <- colnames(pbmc)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCA_common_')
pbmc[["common_factor"]] <- zz

png("../../out/Writeup10_pbmc_10x_dcca_common_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Common UMAP (D-CCA, Meta-cells)")
graphics.off()

###

tmp <- RSpectra::svds(res$common_mat_1 + res$distinct_mat_1, k = K)
tmp <- multiomicCCA:::.mult_mat_vec(tmp$u, tmp$d)
rownames(tmp) <- colnames(pbmc)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCA_common_')
pbmc[["common_factor"]] <- zz

png("../../out/Writeup10_pbmc_10x_dcca_rna_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA-only UMAP (D-CCA, Meta-cells)")
graphics.off()

###

tmp <- RSpectra::svds(res$common_mat_2 + res$distinct_mat_2, k = K)
tmp <- multiomicCCA:::.mult_mat_vec(tmp$u, tmp$d)
rownames(tmp) <- colnames(pbmc)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCA_common_')
pbmc[["common_factor"]] <- zz

png("../../out/Writeup10_pbmc_10x_dcca_atac_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("ATAC-only UMAP (D-CCA, Meta-cells)")
graphics.off()

####

tmp <- RSpectra::svds(res$distinct_mat_1, k = K)
tmp <- multiomicCCA:::.mult_mat_vec(tmp$u, tmp$d)
rownames(tmp) <- colnames(pbmc)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCA_common_')
pbmc[["common_factor"]] <- zz

png("../../out/Writeup10_pbmc_10x_dcca_rna_distinct_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA distinct UMAP (D-CCA, Meta-cells)")
graphics.off()

####

tmp <- RSpectra::svds(res$distinct_mat_2, k = K)
tmp <- multiomicCCA:::.mult_mat_vec(tmp$u, tmp$d)
rownames(tmp) <- colnames(pbmc)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCA_common_')
pbmc[["common_factor"]] <- zz

png("../../out/Writeup10_pbmc_10x_dcca_atac_distinct_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("ATAC distinct UMAP (D-CCA, Meta-cells)")
graphics.off()

