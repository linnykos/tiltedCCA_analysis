rm(list=ls())
library(Seurat)
load("../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_dcca4_zongming.RData")
dcca_res$tilt_perc
dcca_res <- dcca_res2

png("../../../../out/figures/Writeup14h/Writeup14h_citeseq_pbmc224_dcca4_dcca_scores.png",
    height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
multiomicCCA::plot_scores_heatmap.dcca(dcca_res,
                                       membership_vec = as.factor(pbmc$celltype.l2),
                                       log_scale = T, scaling_power = 4)
graphics.off()

png("../../../../out/figures/Writeup14h/Writeup14h_citeseq_pbmc224_dcca4_cca_scores.png",
    height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
multiomicCCA::plot_scores_heatmap.list(list(dcca_res$score_1, dcca_res$score_2),
                                       main_vec = c("CCA scores 1", "CCA scores 2"),
                                       membership_vec = as.factor(pbmc$celltype.l2),
                                       log_scale = T, scaling_power = 4)
graphics.off()

png("../../../../out/figures/Writeup14h/Writeup14h_citeseq_pbmc224_dcca4_svd.png",
    height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
multiomicCCA::plot_scores_heatmap.list(list(multiomicCCA:::.mult_mat_vec(dcca_res$svd_1$u, dcca_res$svd_1$d), 
                                            multiomicCCA:::.mult_mat_vec(dcca_res$svd_2$u, dcca_res$svd_2$d)),
                                       main_vec = c("SVD 1", "SVD 2"),
                                       membership_vec = as.factor(pbmc$celltype.l2),
                                       log_scale = T, scaling_power = 4)
graphics.off()

multiomicCCA:::.l2norm(dcca_res$common_score[,1])
multiomicCCA:::.l2norm(dcca_res$distinct_score_1[,1])
multiomicCCA:::.l2norm(dcca_res$distinct_score_2[,1])

############################

residual_mat_1 <- sapply(1:ncol(dcca_res$distinct_score_1), function(j){
  df <- cbind(dcca_res$distinct_score_1[,j], dcca_res$common_score)
  df <- as.data.frame(df)
  colnames(df) <- c("y", paste0("x", 1:ncol(dcca_res$common_score)))
  lm_res <- stats::lm(y ~ . , data = df)
  stats::residuals(lm_res)
})
alignment_1 <- 1 - sum((residual_mat_1)^2)/sum((dcca_res$distinct_score_1)^2)

residual_mat_2 <- sapply(1:ncol(dcca_res$distinct_score_2), function(j){
  df <- cbind(dcca_res$distinct_score_2[,j], dcca_res$common_score)
  df <- as.data.frame(df)
  colnames(df) <- c("y", paste0("x", 1:ncol(dcca_res$common_score)))
  lm_res <- stats::lm(y ~ . , data = df)
  stats::residuals(lm_res)
})
alignment_2 <- 1 - sum((residual_mat_2)^2)/sum((dcca_res$distinct_score_2)^2)

#######################################
dcca_res <- dcca_res
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res)

n <- nrow(pbmc@meta.data)
svd_1 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_1, 
                                       K = rank_1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_2, K = rank_2, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*sqrt(n))
l2_vec <- apply(dimred_1, 1, function(x){multiomicCCA:::.l2norm(x)})
dimred_1 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, dimred_1)
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*sqrt(n))
l2_vec <- apply(dimred_2, 1, function(x){multiomicCCA:::.l2norm(x)})
dimred_2 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, dimred_2)

set.seed(10)
common_umap <- Seurat::RunUMAP(cbind(dimred_1, dimred_2), 
                               metric = "euclidean",
                               reduction.key = "umapCommon_")
rownames(common_umap@cell.embeddings) <- rownames(pbmc@meta.data)
pbmc[["dcca_common"]] <- Seurat::CreateDimReducObject(common_umap@cell.embeddings, 
                                                      assay = "SCT")

#######################################

svd_1 <- multiomicCCA:::.svd_truncated(dcca_decomp$distinct_mat_1, 
                                       K = rank_1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)

set.seed(10)
distinct1_umap <- Seurat::RunUMAP(dimred_1, 
                                  metric = "cosine",
                                  reduction.key = "umapDistinct1_")
rownames(distinct1_umap@cell.embeddings) <- rownames(pbmc@meta.data)
pbmc[["dcca_distinct1"]] <- Seurat::CreateDimReducObject(distinct1_umap@cell.embeddings, 
                                                         assay = "SCT")

#######################################

svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$distinct_mat_2, 
                                       K = rank_2, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d)

set.seed(10)
distinct2_umap <- Seurat::RunUMAP(dimred_2, 
                                  metric = "cosine",
                                  reduction.key = "umapDistinct2_")
rownames(distinct2_umap@cell.embeddings) <- rownames(pbmc@meta.data)
pbmc[["dcca_distinct2"]] <- Seurat::CreateDimReducObject(distinct2_umap@cell.embeddings,
                                                         assay = "ADT")

#######################################

anchor_name <- "rna.umap"
other_names <- c("dcca_common", "dcca_distinct1", "dcca_distinct2")

for(umap_name in other_names){
  print(umap_name)
  u_mat1 <- pbmc[[anchor_name]]@cell.embeddings
  u_mat2 <- pbmc[[umap_name]]@cell.embeddings
  tmp <- svd(t(u_mat1) %*% u_mat2)
  rotation_mat <- tmp$u %*% t(tmp$v)
  tmp <- u_mat2 %*% t(rotation_mat)
  rownames(tmp) <- rownames(pbmc@meta.data)
  colnames(tmp) <- colnames(pbmc[[umap_name]]@cell.embeddings)
  pbmc[[umap_name]]@cell.embeddings <- tmp
}

# plot according to clones
reduction_vec <- c(anchor_name, other_names)
group_vec <- c("celltype.l2")
main_vec <- c("(RNA)", 
              "(D-CCA, Common)", 
              paste0("(D-CCA, Distinct 1, Alignment: ", round(alignment_1, 2),")"), 
              paste0("(D-CCA, Distinct 2, Alignment: ", round(alignment_2, 2),")"))
file_vec <- c("rna", "dcca-common", "dcca-distinct1", "dcca-distinct2")

for(i in 1:length(reduction_vec)){
  for(j in 1:length(group_vec)){
    plot1 <- Seurat::DimPlot(pbmc, reduction = reduction_vec[i],
                             group.by = group_vec[j], label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\n", main_vec[i]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14h/Writeup14h_citeseq_pbmc224_dcca4_", file_vec[i], ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}
