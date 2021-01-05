rm(list=ls()); set.seed(10); # gcinfo(TRUE)

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
zz <- extract_embedding(res, common_1 = T, common_2 = T, distinct_1 = F, distinct_2 = F, 
                        only_embedding = F)
pbmc[["common_factor"]] <- zz
png("../../out/figures/Writeup10/Writeup10_pbmc_10x_dcca_common_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Common view (D-CCA, Meta-cells)")
graphics.off()

set.seed(10)
zz <- extract_embedding(res, common_1 = F, common_2 = F, distinct_1 = T, distinct_2 = T, 
                        only_embedding = F)
pbmc[["common_factor"]] <- zz
png("../../out/figures/Writeup10/Writeup10_pbmc_10x_dcca_distinct_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Distinct view (D-CCA, Meta-cells)")
graphics.off()

set.seed(10)
zz <- extract_embedding(res, common_1 = T, common_2 = T, distinct_1 = T, distinct_2 = T, only_embedding = F)
pbmc[["common_factor"]] <- zz
png("../../out/figures/Writeup10/Writeup10_pbmc_10x_dcca_everything_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Everything view (D-CCA, Meta-cells)")
graphics.off()

##################################

set.seed(10)
zz <- extract_embedding(res, common_1 = T, common_2 = F, distinct_1 = F, distinct_2 = F, 
                        only_embedding = F)
pbmc[["common_factor"]] <- zz
png("../../out/figures/Writeup10/Writeup10_pbmc_10x_dcca_rna_common_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA common view (D-CCA, Meta-cells)")
graphics.off()

set.seed(10)
zz <- extract_embedding(res, common_1 = F, common_2 = F, distinct_1 = T, distinct_2 = F, 
                        only_embedding = F)
pbmc[["common_factor"]] <- zz
png("../../out/figures/Writeup10/Writeup10_pbmc_10x_dcca_rna_common_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA distinct view (D-CCA, Meta-cells)")
graphics.off()

set.seed(10)
zz <- extract_embedding(res, common_1 = T, common_2 = F, distinct_1 = T, distinct_2 = F, 
                        only_embedding = F)
pbmc[["common_factor"]] <- zz
png("../../out/figures/Writeup10/Writeup10_pbmc_10x_dcca_rna_only_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA view (D-CCA, Meta-cells)")
graphics.off()

###

set.seed(10)
zz <- extract_embedding(res, common_1 = F, common_2 = T, distinct_1 = F, distinct_2 = F, 
                        only_embedding = F)
pbmc[["common_factor"]] <- zz
png("../../out/figures/Writeup10/Writeup10_pbmc_10x_dcca_atac_common_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("ATAC common view (D-CCA, Meta-cells)")
graphics.off()

set.seed(10)
zz <- extract_embedding(res, common_1 = F, common_2 = F, distinct_1 = F, distinct_2 = T, 
                        only_embedding = F)
pbmc[["common_factor"]] <- zz
png("../../out/figures/Writeup10/Writeup10_pbmc_10x_dcca_atac_common_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("ATAC distinct view (D-CCA, Meta-cells)")
graphics.off()

set.seed(10)
zz <- extract_embedding(res, common_1 = F, common_2 = T, distinct_1 = F, distinct_2 = T, 
                        only_embedding = F)
pbmc[["common_factor"]] <- zz
png("../../out/figures/Writeup10/Writeup10_pbmc_10x_dcca_atac_only_umap_meta.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'common_factor', group.by = 'predicted.id', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("ATAC view (D-CCA, Meta-cells)")
graphics.off()

