rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)
library(multiomicCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat_1 <- t(bm[["RNA"]]@scale.data)
mat_2 <- t(bm[["ADT"]]@scale.data)

set.seed(10)
rank_1 <- 30; rank_2 <- 18
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, 
                                      rank_1 = rank_1, rank_2 = rank_2, 
                                      apply_shrinkage = F, verbose = T) 
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, rank_c = min(rank_1, rank_2), verbose = T)
svd_list <- multiomicCCA::extract_svd_embedding(dcca_decomp)

########################

set.seed(10)
zz <- multiomicCCA::extract_umap_embedding(svd_list, common_1 = T, common_2 = T, distinct_1 = F, distinct_2 = F, 
                        only_embedding = F, vis_param = dcca_decomp$vis_param)
bm[["common_factor"]] <- zz
png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_common_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm, reduction = 'common_factor', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Common view (25, D-CCA)")
graphics.off()

set.seed(10)
zz <- multiomicCCA::extract_umap_embedding(svd_list, common_1 = F, common_2 = F, distinct_1 = T, distinct_2 = T, 
                                      only_embedding = F, vis_param = dcca_decomp$vis_param)
bm[["common_factor"]] <- zz
png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_distinct_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm, reduction = 'common_factor', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Distinct view (25, D-CCA)")
graphics.off()

set.seed(10)
zz <- multiomicCCA::extract_umap_embedding(svd_list, common_1 = T, common_2 = T, distinct_1 = T, distinct_2 = T, 
                                      only_embedding = F, vis_param = dcca_decomp$vis_param)
bm[["common_factor"]] <- zz
png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_everything_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm, reduction = 'common_factor', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Everything view (25, D-CCA)")
graphics.off()

##################

set.seed(10)
zz <- multiomicCCA::extract_umap_embedding(svd_list, common_1 = T, common_2 = F, distinct_1 = F, distinct_2 = F, 
                                      only_embedding = F, vis_param = dcca_decomp$vis_param)
bm[["common_factor"]] <- zz
png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_rna_common_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm, reduction = 'common_factor', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA common view (25, D-CCA)")
graphics.off()

set.seed(10)
zz <- multiomicCCA::extract_umap_embedding(svd_list, common_1 = F, common_2 = F, distinct_1 = T, distinct_2 = F, 
                                      only_embedding = F, vis_param = dcca_decomp$vis_param)
bm[["common_factor"]] <- zz
png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_rna_distinct_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm, reduction = 'common_factor', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA distinct view (25, D-CCA)")
graphics.off()

set.seed(10)
zz <- multiomicCCA::extract_umap_embedding(svd_list, common_1 = T, common_2 = F, distinct_1 = T, distinct_2 = F, 
                                      only_embedding = F, vis_param = dcca_decomp$vis_param)
bm[["common_factor"]] <- zz
png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_rna_everything_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm, reduction = 'common_factor', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA everything view (25, D-CCA)")
graphics.off()

##################

set.seed(10)
zz <- multiomicCCA::extract_umap_embedding(svd_list, common_1 = F, common_2 = T, distinct_1 = F, distinct_2 = F, 
                                      only_embedding = F, vis_param = dcca_decomp$vis_param)
bm[["common_factor"]] <- zz
png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_protein_common_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm, reduction = 'common_factor', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Protein common view (25, D-CCA)")
graphics.off()

set.seed(10)
zz <- multiomicCCA::extract_umap_embedding(svd_list, common_1 = F, common_2 = F, distinct_1 = F, distinct_2 = T, 
                                      only_embedding = F, vis_param = dcca_decomp$vis_param)
bm[["common_factor"]] <- zz
png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_protein_distinct_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm, reduction = 'common_factor', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Protein distinct view (25, D-CCA)")
graphics.off()

set.seed(10)
zz <- multiomicCCA::extract_umap_embedding(svd_list, common_1 = F, common_2 = T, distinct_1 = F, distinct_2 = T, 
                                      only_embedding = F, vis_param = dcca_decomp$vis_param)
bm[["common_factor"]] <- zz
png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_protein_everything_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm, reduction = 'common_factor', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Protein everything view (25, D-CCA)")
graphics.off()

#############################

set.seed(10)
tmp <- multiomicCCA::dcca_information_weight(mat = dcca_decomp$common_mat_1+dcca_decomp$distinct_mat_1,
                               common_mat = dcca_decomp$common_mat_1, verbose = T, iter_max = 20)
weight_vec1 <- tmp$alpha_vec
tmp <- multiomicCCA::dcca_information_weight(mat = dcca_decomp$common_mat_2+dcca_decomp$distinct_mat_2,
                               common_mat = dcca_decomp$common_mat_2, verbose = T, iter_max = 20)
weight_vec2 <- tmp$alpha_vec
weight_mat <- data.frame(RNA = 1-weight_vec1, Protein = 1-weight_vec2, cell_type = bm@meta.data$celltype.l2)
png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_information_rank1.png", height = 1500, width = 2500, units = "px", res = 300)
par(mfrow = c(1,1), mar = c(6, 4, 3.5, 0.5))
multiomicCCA::information_plot(weight_mat,  main = "D-CCA's weights (25, Rank-1)")
graphics.off()

mean(weight_vec1); mean(weight_vec2)
mean(c(weight_vec1, weight_vec2))

custom_explained_variance <- function(dcca_decomp, weight_func = function(x){x}){
  stopifnot(class(dcca_decomp) == "dcca_decomp")
  
  n <- nrow(dcca_decomp$common_mat_1)
  p1 <- ncol(dcca_decomp$common_mat_1); p2 <- ncol(dcca_decomp$common_mat_2)
  
  cell_weight_vec1 <- sapply(1:n, function(i){
    val1 <- weight_func(multiomicCCA:::.l2norm(dcca_decomp$common_mat_1[i,]))
    val2 <- weight_func(multiomicCCA:::.l2norm(dcca_decomp$distinct_mat_1[i,]))
    val2/(val1+val2)
  })
  
  cell_weight_vec2 <- sapply(1:n, function(i){
    val1 <- weight_func(multiomicCCA:::.l2norm(dcca_decomp$common_mat_2[i,]))
    val2 <- weight_func(multiomicCCA:::.l2norm(dcca_decomp$distinct_mat_2[i,]))
    val2/(val1+val2)
  })
  
  list(cell_weight_vec1 = cell_weight_vec1,
       cell_weight_vec2 = cell_weight_vec2)
}

custom_explained_variance2 <- function(dcca_decomp, membership_vec){
  stopifnot(class(dcca_decomp) == "dcca_decomp")
  
  n <- nrow(dcca_decomp$common_mat_1)
  mat_1 <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
  weight_vec_1 <- sapply(sort(unique(membership_vec)), function(i){
    idx <- which(membership_vec == i)
    val <- multiomicCCA:::.l2norm(dcca_decomp$distinct_mat_1[idx,])^2/multiomicCCA:::.l2norm(mat_1[idx,])^2
    max(min(val, 1),0)
  })
  
  mat_2 <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2
  weight_vec_2 <- sapply(sort(unique(membership_vec)), function(i){
    idx <- which(membership_vec == i)
    val <- multiomicCCA:::.l2norm(dcca_decomp$distinct_mat_2[idx,])^2/multiomicCCA:::.l2norm(mat_2[idx,])^2
    max(min(val, 1),0)
  })
  
  list(weight_vec_1 = weight_vec_1,
       weight_vec_2 = weight_vec_2)
}

tmp <- custom_explained_variance2(dcca_decomp, bm@meta.data$celltype.l2)
weight_mat <- data.frame(RNA = tmp$weight_vec_1, Protein = tmp$weight_vec_2, cell_type = sort(unique(bm@meta.data$celltype.l2)))
png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_information_R2.png", height = 1500, width = 2500, units = "px", res = 300)
par(mfrow = c(1,1), mar = c(6, 4, 3.5, 0.5))
multiomicCCA::information_plot(weight_mat, main = "D-CCA's weights (25, R-Squared)", plot_individual = F)
graphics.off()

tmp <- custom_explained_variance(dcca_decomp)
weight_mat <- data.frame(RNA = tmp$cell_weight_vec1, Protein = tmp$cell_weight_vec2, cell_type = bm@meta.data$celltype.l2)
png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_information_ratio.png", height = 1500, width = 2500, units = "px", res = 300)
par(mfrow = c(1,1), mar = c(6, 4, 3.5, 0.5))
multiomicCCA::information_plot(weight_mat, main = "D-CCA's weights (25, L2 ratio)")
graphics.off()

tmp <- custom_explained_variance(dcca_decomp, weight_func = function(x){x^2})
weight_mat <- data.frame(RNA = tmp$cell_weight_vec1, Protein = tmp$cell_weight_vec2, cell_type = bm@meta.data$celltype.l2)
png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_information_ratio2.png", height = 1500, width = 2500, units = "px", res = 300)
par(mfrow = c(1,1), mar = c(6, 4, 3.5, 0.5))
multiomicCCA::information_plot(weight_mat, main = "D-CCA's weights (25, Squared L2 ratio)")
graphics.off()


###############################
tmp1 <- sapply(1:ncol(dcca_decomp$common_score), function(j){
  multiomicCCA:::.l2norm(dcca_decomp$common_score[,j])/(multiomicCCA:::.l2norm(dcca_decomp$common_score[,j]) + multiomicCCA:::.l2norm(dcca_decomp$distinct_score_1[,j]))
})
tmp4 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1, rank_1)
tmp5 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2, rank_2)

png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_dcca_values.png", height = 800, width = 2000, units = "px", res = 300)
par(mfrow = c(1,3))
# cca values
plot(dcca_res$cca_obj, main = "CCA objective", xlab = "Canonical dimension", 
     ylab = "Canonical correlation", pch = 16)

# amount of variability, common vs distincts
plot(NA, xlim = c(1,length(tmp1)), ylim = c(0,1), main = "Ratio of canonical\nscores' signal",
     xlab = "Canonical dimension", ylab = "Ratio of signals (common/total)")
points(tmp1, pch = 16)
for(i in 1:length(tmp1)){
  lines(rep(i,2), c(0,tmp1[i]), col = 1)
  lines(rep(i,2), c(1,tmp1[i]), col = 2)
}
legend("topright", c("Common", "Distinct"), fill= c(1,2))

# amount of variability, rna vs atac
ylim <- range(c(tmp4$d, tmp5$d))
plot(NA, xlim = c(1,max(c(rank_1,rank_2))), ylim = ylim, main = "Singular value of modalities",
     xlab = "Singular dimension", ylab = "Singular value")
points(tmp4$d, pch = 16)
points(tmp5$d, col = 2, pch = 16)
legend("topright", c("RNA", "Protein"), bty="n", fill=c(1,2))
graphics.off()



