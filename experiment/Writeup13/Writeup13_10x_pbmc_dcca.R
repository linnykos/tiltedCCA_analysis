rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../../../out/Writeup10/Writeup10_10x_pbmc_preprocess4.RData")
date_of_run <- Sys.time(); session_info <- sessionInfo()

mat_1 <- t(pbmc[["SCT"]]@scale.data)
mat_2 <- t(pbmc[["ATAC"]]@scale.data)

head(rownames(mat_1)); head(colnames(mat_1))
head(rownames(mat_2)); head(colnames(mat_2))
dim(mat_1); dim(mat_2)
metadata <- pbmc@meta.data
table(metadata$predicted.id)
rm(list = "pbmc"); gc(T)

set.seed(10)
rank_1 <- 17; rank_2 <- 13
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, rank_1 = rank_1, rank_2 = rank_2, 
                                      meta_clustering = NA, num_neigh = 100, 
                                      apply_shrinkage = F, fix_distinct_perc = F, 
                                      verbose = T) 
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, rank_c = min(rank_1, rank_2), 
                                                verbose = T)
dcca_decomp$cca_obj
dcca_decomp$distinct_perc_2
rm(list = c("mat_1", "mat_2")); gc(T)

save.image("../../../../out/Writeup13/Writeup13_10x_pbmc_dcca.RData")

##################


membership_vec <- as.factor(metadata$predicted.id)
png(paste0("../../../../out/figures/Writeup13/Writeup13_10x_pbmc_dcca_scores.png"), height = 2000, width = 2400, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
multiomicCCA::plot_scores(dcca_decomp, membership_vec, decomposition = F)
graphics.off()

png(paste0("../../../../out/figures/Writeup13/Writeup13_10x_pbmc_dcca_scores_decomp.png"), height = 2000, width = 3600, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
multiomicCCA::plot_scores(dcca_decomp, membership_vec, decomposition = T)
graphics.off()

png(paste0("../../../../out/figures/Writeup13/Writeup13_10x_pbmc_dcca_score_heatmap.png"), height = 1200, width = 2400, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
multiomicCCA::plot_scores_heatmap(dcca_decomp, membership_vec, log_scale = T, scaling_power = 1.5)
graphics.off()

png(paste0("../../../../out/figures/Writeup13/Writeup13_10x_pbmc_dcca_summary.png"), height = 1200, width = 1200, units = "px", res = 300)
par(mar = c(4,4,4,4))
multiomicCCA::plot_summary(dcca_decomp, main = "Summary of D-CCA")
graphics.off()


##################

load("../../../../out/Writeup10/Writeup10_10x_pbmc_preprocess4.RData")
Seurat::DefaultAssay(pbmc) <- "SCT"
pbmc[["ATAC"]] <- NULL; pbmc[["wsnn"]] <- NULL;
pbmc[["pca"]] <- NULL; pbmc[["umap.rna"]] <- NULL; pbmc[["umap.atac"]] <- NULL
pbmc[["RNA"]] <- NULL; pbmc[["wknn"]] <- NULL; pbmc[["wnn.umap"]] <- NULL

membership_vec <- as.factor(metadata$predicted.id)
main_vec <- c("common", "distinct", "everything")
set.seed(10)
zz <- multiomicCCA::plot_embeddings(dcca_decomp, membership_vec, data_1 = T, data_2 = F, 
                                    add_noise = T, pca = F, only_embedding = T, verbose = T)

for(i in 1:3){
  pbmc[["asdf"]] <- Seurat::CreateDimReducObject(embedding = zz[[i]], key = "UMAP", assay = "SCT")
  png(paste0("../../../../out/figures/Writeup13/Writeup13_10x_pbmc_dcca_rna_", main_vec[i], "_umap.png"), height = 1500, width = 1500, units = "px", res = 300)
  plot1 <- Seurat::DimPlot(pbmc, reduction = 'asdf', group.by = 'predicted.id', label = TRUE, 
                           repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
  plot1 + ggplot2::ggtitle(paste0("RNA ", main_vec[i], "  view (10x, D-CCA)"))
  graphics.off()
}

#####

membership_vec <- as.factor(metadata$predicted.id)
main_vec <- c("common", "distinct", "everything")
set.seed(10)
zz <- multiomicCCA::plot_embeddings(dcca_decomp, membership_vec, data_1 = F, data_2 = T, 
                                    add_noise = T, pca = F, only_embedding = T, verbose = T)

for(i in 1:3){
  pbmc[["asdf"]] <- Seurat::CreateDimReducObject(embedding = zz[[i]], key = "UMAP", assay = "SCT")
  png(paste0("../../../../out/figures/Writeup13/Writeup13_10x_pbmc_dcca_atac_", main_vec[i], "_umap.png"), height = 1500, width = 1500, units = "px", res = 300)
  plot1 <- Seurat::DimPlot(pbmc, reduction = 'asdf', group.by = 'predicted.id', label = TRUE, 
                           repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
  plot1 + ggplot2::ggtitle(paste0("ATAC ", main_vec[i], "  view (10x, D-CCA)"))
  graphics.off()
}

#####

membership_vec <- as.factor(metadata$predicted.id)
main_vec <- c("common", "distinct", "everything")
set.seed(10)
zz <- multiomicCCA::plot_embeddings(dcca_decomp, membership_vec, data_1 = T, data_2 = T, 
                                    add_noise = T, pca = F, only_embedding = T, verbose = T)

for(i in 1:3){
  pbmc[["asdf"]] <- Seurat::CreateDimReducObject(embedding = zz[[i]], key = "UMAP", assay = "SCT")
  png(paste0("../../../../out/figures/Writeup13/Writeup13_10x_pbmc_dcca_both_", main_vec[i], "_umap.png"), height = 1500, width = 1500, units = "px", res = 300)
  plot1 <- Seurat::DimPlot(pbmc, reduction = 'asdf', group.by = 'predicted.id', label = TRUE, 
                           repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
  plot1 + ggplot2::ggtitle(paste0("Both ", main_vec[i], "  view (10x, D-CCA)"))
  graphics.off()
}

#######################

membership_vec <- as.factor(metadata$predicted.id)
celltype_idx <- as.numeric(which(table(membership_vec) <= 10))
celltype_rm <- names(which(table(membership_vec) <= 10))
idx <- which(membership_vec %in% celltype_rm)

membership_vec2 <- as.factor(as.character(membership_vec[-idx]))
set.seed(10)
clisi_1 <- multiomicCCA::clisi_information(dcca_decomp$common_mat_1[-idx,,drop=F], dcca_decomp$distinct_mat_1[-idx,,drop=F],
                                           membership_vec = membership_vec2, 
                                           rank_c = min(rank_1, rank_2), rank_d = rank_1, 
                                           nn = 100, radius_quantile = 0.9,
                                           max_subsample_clisi = 1000,
                                           verbose = T)
set.seed(10)
clisi_2 <- multiomicCCA::clisi_information(dcca_decomp$common_mat_2[-idx,,drop=F], dcca_decomp$distinct_mat_2[-idx,,drop=F],
                                           membership_vec = membership_vec2, 
                                           rank_c = min(rank_1, rank_2), rank_d = rank_2, 
                                           nn = 100, radius_quantile = 0.9,
                                           max_subsample_clisi = 1000,
                                           verbose = T)

png(paste0("../../../../out/figures/Writeup13/Writeup13_10x_pbmc_clisi.png"), height = 1500, width = 3000, units = "px", res = 300)
multiomicCCA::plot_clisi(clisi_1, clisi_2, col_vec = scales::hue_pal()(length(levels(membership_vec)))[-celltype_idx],
                         main = "cLISI (10x, D-CCA)")
graphics.off()

png(paste0("../../../../out/figures/Writeup13/Writeup13_10x_pbmc_clisi_legend.png"), height = 1500, width = 1000, units = "px", res = 300)
multiomicCCA:::plot_clisi_legend(clisi_1, col_vec = scales::hue_pal()(length(levels(membership_vec)))[-celltype_idx])
graphics.off()

cbind(as.character(clisi_1$common_clisi$membership_info[,1]), round(clisi_1$common_clisi$membership_info[,4],2), round(clisi_1$distinct$membership_info[,4],2))
cbind(as.character(clisi_2$common_clisi$membership_info[,1]), round(clisi_2$common_clisi$membership_info[,4],2), round(clisi_2$distinct$membership_info[,4],2))

save.image("../../../../out/Writeup13/Writeup13_10x_pbmc_dcca.RData")

###############################3
################################

load("../../../../out/Writeup10/Writeup10_10x_pbmc_preprocess4.RData")
Seurat::DefaultAssay(pbmc) <- "SCT"
pbmc[["ATAC"]] <- NULL; pbmc[["wsnn"]] <- NULL;
pbmc[["pca"]] <- NULL; pbmc[["umap.rna"]] <- NULL; pbmc[["umap.atac"]] <- NULL
pbmc[["RNA"]] <- NULL; pbmc[["wknn"]] <- NULL; pbmc[["wnn.umap"]] <- NULL

membership_vec <- as.factor(metadata$predicted.id)
main_vec <- c("common", "distinct", "everything")
set.seed(10)
zz <- multiomicCCA::plot_embeddings(dcca_decomp, membership_vec, data_1 = T, data_2 = F, 
                                    add_noise = F, pca = F, only_embedding = T, verbose = T)

for(i in 1:3){
  pbmc[["asdf"]] <- Seurat::CreateDimReducObject(embedding = zz[[i]], key = "UMAP", assay = "SCT")
  png(paste0("../../../../out/figures/Writeup13/Writeup13_10x_pbmc_dcca_rna_", main_vec[i], "_umap_nonoise.png"), height = 1500, width = 1500, units = "px", res = 300)
  plot1 <- Seurat::DimPlot(pbmc, reduction = 'asdf', group.by = 'predicted.id', label = TRUE, 
                           repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
  plot1 + ggplot2::ggtitle(paste0("RNA ", main_vec[i], "  view (10x, D-CCA)"))
  graphics.off()
}

membership_vec <- as.factor(metadata$predicted.id)
main_vec <- c("common", "distinct", "everything")
set.seed(10)
zz <- multiomicCCA::plot_embeddings(dcca_decomp, membership_vec, data_1 = F, data_2 = T, 
                                    add_noise = F, pca = F, only_embedding = T, verbose = T)

for(i in 1:3){
  pbmc[["asdf"]] <- Seurat::CreateDimReducObject(embedding = zz[[i]], key = "UMAP", assay = "SCT")
  png(paste0("../../../../out/figures/Writeup13/Writeup13_10x_pbmc_dcca_atac_", main_vec[i], "_umap_nonoise.png"), height = 1500, width = 1500, units = "px", res = 300)
  plot1 <- Seurat::DimPlot(pbmc, reduction = 'asdf', group.by = 'predicted.id', label = TRUE, 
                           repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
  plot1 + ggplot2::ggtitle(paste0("ATAC ", main_vec[i], "  view (10x, D-CCA)"))
  graphics.off()
}
