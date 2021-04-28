rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_pbmc228_preprocessed2.RData")

library(Seurat)
library(multiomicCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat_1 <- t(pbmc[["SCT"]]@scale.data)
mat_2 <- t(pbmc[["ADT"]]@scale.data)
dim(mat_1); dim(mat_2)
metadata <- pbmc@meta.data
rm(list = "pbmc"); gc(T)

set.seed(10)
rank_1 <- 40; rank_2 <- 50
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, rank_1 = rank_1, rank_2 = rank_2, 
                                      meta_clustering = NA, num_neigh = 100, cell_max = 15000,
                                      apply_shrinkage = F, fix_distinct_perc = F, 
                                      verbose = T) 
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, rank_c = min(rank_1, rank_2), 
                                                verbose = T)
dcca_decomp$cca_obj
dcca_decomp$distinct_perc_2
rm(list = c("mat_1", "mat_2", "zz", "plot1")); gc(T)

save.image("../../../../out/Writeup13/Writeup13_citeseq_pbmc224_dcca.RData")

###########

membership_vec <- as.factor(metadata$celltype.l2)
png(paste0("../../../../out/figures/Writeup13/Writeup13_citeseq_pbmc224_dcca_scores.png"), height = 2000, width = 2400, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
multiomicCCA::plot_scores(dcca_decomp, membership_vec, decomposition = F)
graphics.off()

png(paste0("../../../../out/figures/Writeup13/Writeup13_citeseq_pbmc224_dcca_scores_decomp.png"), height = 2000, width = 3600, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
multiomicCCA::plot_scores(dcca_decomp, membership_vec, decomposition = T)
graphics.off()

png(paste0("../../../../out/figures/Writeup13/Writeup13_citeseq_pbmc224_dcca_score_heatmap.png"), height = 1200, width = 2400, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
multiomicCCA::plot_scores_heatmap(dcca_decomp, membership_vec, log_scale = T, scaling_power = 2)
graphics.off()

png(paste0("../../../../out/figures/Writeup13/Writeup13_citeseq_pbmc224_dcca_summary.png"), height = 1200, width = 1200, units = "px", res = 300)
par(mar = c(4,4,4,4))
multiomicCCA::plot_summary(dcca_decomp, main = "Summary of D-CCA")
graphics.off()


##########

load("../../../../out/Writeup11/Writeup11_citeseq_pbmc228_preprocessed2.RData")
Seurat::DefaultAssay(pbmc) <- "ADT"
pbmc[["SCT"]] <- NULL; pbmc[["wsnn"]] <- NULL;
pbmc[["apca"]] <- NULL; pbmc[["aumap"]] <- NULL; pbmc[["pca"]] <- NULL; pbmc[["spca"]] <- NULL; 
pbmc[["umap"]] <- NULL; pbmc[["rna.umap"]] <- NULL; pbmc[["adt.umap"]] <- NULL; pbmc[["wnn.umap"]] <- NULL; 

membership_vec <- as.factor(metadata$celltype.l2) 
main_vec <- c("common", "distinct", "everything")
set.seed(10)
zz <- multiomicCCA::plot_embeddings(dcca_decomp, membership_vec, data_1 = T, data_2 = F, 
                                    add_noise = T, pca = F, only_embedding = T, verbose = T)
for(i in 1:3){
  pbmc[["asdf"]] <- Seurat::CreateDimReducObject(embedding = zz[[i]], key = "UMAP", assay = "ADT")
  png(paste0("../../../../out/figures/Writeup13/Writeup13_citeseq_pbmc224_dcca_rna_", main_vec[i], "_umap.png"), height = 1500, width = 1500, units = "px", res = 300)
  plot1 <- Seurat::DimPlot(pbmc, reduction = 'asdf', group.by = 'celltype.l2', label = TRUE, 
                           repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
  plot1 + ggplot2::ggtitle(paste0("RNA ", main_vec[i], "  view (224, D-CCA)"))
  graphics.off()
}

#####

membership_vec <- as.factor(metadata$celltype.l2) 
main_vec <- c("common", "distinct", "everything")
set.seed(10)
zz <- multiomicCCA::plot_embeddings(dcca_decomp, membership_vec, data_1 = F, data_2 = T, 
                                    add_noise = T, pca = F, only_embedding = T, verbose = T)
for(i in 1:3){
  pbmc[["asdf"]] <- Seurat::CreateDimReducObject(embedding = zz[[i]], key = "UMAP", assay = "ADT")
  png(paste0("../../../../out/figures/Writeup13/Writeup13_citeseq_pbmc224_dcca_protein_", main_vec[i], "_umap.png"), height = 1500, width = 1500, units = "px", res = 300)
  plot1 <- Seurat::DimPlot(pbmc, reduction = 'asdf', group.by = 'celltype.l2', label = TRUE, 
                           repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
  plot1 + ggplot2::ggtitle(paste0("Protein ", main_vec[i], "  view (224, D-CCA)"))
  graphics.off()
}

#####

membership_vec <- as.factor(metadata$celltype.l2) 
main_vec <- c("common", "distinct", "everything")
set.seed(10)
zz <- multiomicCCA::plot_embeddings(dcca_decomp, membership_vec, data_1 = T, data_2 = T, 
                                    add_noise = T, pca = F, only_embedding = T, verbose = T)
for(i in 1:3){
  pbmc[["asdf"]] <- Seurat::CreateDimReducObject(embedding = zz[[i]], key = "UMAP", assay = "ADT")
  png(paste0("../../../../out/figures/Writeup13/Writeup13_citeseq_pbmc224_dcca_both_", main_vec[i], "_umap.png"), height = 1500, width = 1500, units = "px", res = 300)
  plot1 <- Seurat::DimPlot(pbmc, reduction = 'asdf', group.by = 'celltype.l2', label = TRUE, 
                           repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
  plot1 + ggplot2::ggtitle(paste0("Both ", main_vec[i], "  view (224, D-CCA)"))
  graphics.off()
}

#################

membership_vec <- as.factor(metadata$celltype.l2)
set.seed(10)
clisi_1 <- multiomicCCA::clisi_information(dcca_decomp$common_mat_1, dcca_decomp$distinct_mat_1,
                                           membership_vec = membership_vec, 
                                           rank_c = min(rank_1, rank_2), rank_d = rank_1, 
                                           nn = 500, radius_quantile = 0.9,
                                           max_subsample_frnn = 5000, max_subsample_clisi = 1000,
                                           verbose = T)
set.seed(10)
clisi_2 <- multiomicCCA::clisi_information(dcca_decomp$common_mat_2, dcca_decomp$distinct_mat_2,
                                           membership_vec = membership_vec, 
                                           rank_c = min(rank_1, rank_2), rank_d = rank_2, 
                                           nn = 500, radius_quantile = 0.9,
                                           max_subsample_frnn = 5000, max_subsample_clisi = 1000,
                                           verbose = T)

png(paste0("../../../../out/figures/Writeup13/Writeup13_citeseq_pbmc224_clisi.png"), height = 1500, width = 3000, units = "px", res = 300)
multiomicCCA::plot_clisi(clisi_1, clisi_2, main = "cLISI (224, D-CCA)")
graphics.off()

png(paste0("../../../../out/figures/Writeup13/Writeup13_citeseq_pbmc224_clisi_legend.png"), height = 1500, width = 1000, units = "px", res = 300)
multiomicCCA:::plot_clisi_legend(clisi_1)
graphics.off()

cbind(as.character(clisi_1$common_clisi$membership_info[,1]), round(clisi_1$common_clisi$membership_info[,4],2), round(clisi_1$distinct$membership_info[,4],2))
cbind(as.character(clisi_2$common_clisi$membership_info[,1]), round(clisi_2$common_clisi$membership_info[,4],2), round(clisi_2$distinct$membership_info[,4],2))


save.image("../../../../out/Writeup13/Writeup13_citeseq_pbmc224_dcca.RData")

########33333333######################
########33333333######################
########33333333######################

load("../../../../out/Writeup11/Writeup11_citeseq_pbmc228_preprocessed2.RData")
Seurat::DefaultAssay(pbmc) <- "ADT"
pbmc[["SCT"]] <- NULL; pbmc[["wsnn"]] <- NULL;
pbmc[["apca"]] <- NULL; pbmc[["aumap"]] <- NULL; pbmc[["pca"]] <- NULL; pbmc[["spca"]] <- NULL; 
pbmc[["umap"]] <- NULL; pbmc[["rna.umap"]] <- NULL; pbmc[["adt.umap"]] <- NULL; pbmc[["wnn.umap"]] <- NULL; 

membership_vec <- as.factor(metadata$celltype.l2) 
main_vec <- c("common", "distinct", "everything")
set.seed(10)
zz <- multiomicCCA::plot_embeddings(dcca_decomp, membership_vec, data_1 = T, data_2 = F, 
                                    add_noise = F, pca = F, only_embedding = T, verbose = T)



