rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")
date_of_run <- Sys.time(); session_info <- sessionInfo()

# Seurat::DefaultAssay(mbrain) <- "ATAC"
# set.seed(10)
# mbrain <- Seurat::ScaleData(mbrain, features = Seurat::VariableFeatures(object = mbrain))

mat_1 <- t(mbrain[["SCT"]]@scale.data)
Seurat::DefaultAssay(mbrain) <- "ATAC"
mat_2 <- Matrix::t(mbrain[["ATAC"]]@data[Seurat::VariableFeatures(object = mbrain),])

head(rownames(mat_1)); head(colnames(mat_1))
head(rownames(mat_2)); head(colnames(mat_2))
dim(mat_1); dim(mat_2)
metadata <- mbrain@meta.data

# zz <- irlba::irlba(mat_1, nv = 50); zz$d; diff(zz$d)/zz$d[-1]
# zz <- irlba::irlba(mat_2, nv = 50); zz$d; diff(zz$d)/zz$d[-1]

set.seed(10)
rank_1 <- 25; rank_2 <- 25
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, dims_1 = 1:rank_1, dims_2 = 2:rank_2,
                                      center_1 = T, center_2 = T,
                                      meta_clustering = NA, num_neigh = 100, 
                                      apply_shrinkage = F, fix_distinct_perc = F, 
                                      verbose = T) 
dcca_res$cca_obj
dcca_res$distinct_perc_2
dcca_res2 <- dcca_res
class(dcca_res2) <- "dcca_decomp"

# dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res,  verbose = T)
# dcca_decomp$cca_obj
# dcca_decomp$distinct_perc_2
# rm(list = c("mat_1", "mat_2", "mbrain"))
# save.image("../../../../out/Writeup14/Writeup14_10x_mouseembryo_dcca.RData")

##################

membership_vec <- as.factor(metadata$label_Savercat)
png(paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_scores.png"), height = 2000, width = 2400, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
multiomicCCA::plot_scores(dcca_res2, membership_vec, decomposition = F)
graphics.off()

png(paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_scores_decomp.png"), height = 2000, width = 3600, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
multiomicCCA::plot_scores(dcca_res2, membership_vec, decomposition = T)
graphics.off()

png(paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_score_heatmap.png"), height = 1200, width = 2400, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
multiomicCCA::plot_scores_heatmap(dcca_res2, membership_vec, log_scale = T, scaling_power = 1.5)
graphics.off()

png(paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_summary.png"), height = 1200, width = 1200, units = "px", res = 300)
par(mar = c(4,4,4,4))
multiomicCCA::plot_summary(dcca_res2, main = "Summary of D-CCA")
graphics.off()

#####################

Seurat::DefaultAssay(mbrain) <- "SCT"
mbrain[["ATAC"]] <- NULL; mbrain[["wsnn"]] <- NULL;
mbrain[["pca"]] <- NULL; mbrain[["umap"]] <- NULL; mbrain[["umap.atac"]] <- NULL
mbrain[["RNA"]] <- NULL; mbrain[["wknn"]] <- NULL; mbrain[["umap.wnn"]] <- NULL

membership_vec <- as.factor(metadata$label_Savercat)
title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")
set.seed(10)
zz <- multiomicCCA::plot_embeddings(dcca_res2, membership_vec, data_1 = T, data_2 = F, 
                                    add_noise = T, pca = F, only_embedding = T, verbose = T)

for(i in 1:3){
  mbrain[["asdf"]] <- Seurat::CreateDimReducObject(embedding = zz[[i]], key = "UMAP", assay = "SCT")
  plot1 <- Seurat::DimPlot(mbrain, reduction = "asdf", group.by = "label_Savercat", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo (RNA)\n", title_vec[i]) 
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_rna_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 9, height = 5, units = "in")
}

###############

membership_vec <- as.factor(metadata$label_Savercat)
title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")
set.seed(10)
zz <- multiomicCCA::plot_embeddings(dcca_res2, membership_vec, data_1 = F, data_2 = T, 
                                    add_noise = T, pca = F, only_embedding = T, verbose = T)

for(i in 1:3){
  mbrain[["asdf"]] <- Seurat::CreateDimReducObject(embedding = zz[[i]], key = "UMAP", assay = "SCT")
  plot1 <- Seurat::DimPlot(mbrain, reduction = "asdf", group.by = "label_Savercat", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo (ATAC)\n", title_vec[i]) 
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_atac_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 9, height = 5, units = "in")
}

###############

membership_vec <- as.factor(metadata$label_Savercat)
title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")
set.seed(10)
zz <- multiomicCCA::plot_embeddings(dcca_res2, membership_vec, data_1 = T, data_2 = T, 
                                    add_noise = T, pca = F, only_embedding = T, verbose = T)

for(i in 1:3){
  mbrain[["asdf"]] <- Seurat::CreateDimReducObject(embedding = zz[[i]], key = "UMAP", assay = "SCT")
  plot1 <- Seurat::DimPlot(mbrain, reduction = "asdf", group.by = "label_Savercat", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo (ATAC)\n", title_vec[i]) 
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14/Writeup14_10x_mouseembryo_dcca_both_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 9, height = 5, units = "in")
}