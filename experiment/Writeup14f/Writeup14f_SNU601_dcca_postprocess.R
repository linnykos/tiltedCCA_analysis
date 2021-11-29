rm(list=ls())
library(Seurat); library(Signac)
load("../../../../out/Writeup14f/Writeup14f_SNU601_dcca.RData")
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2)

residual_mat_1 <- sapply(1:ncol(dcca_res2$distinct_score_1), function(j){
  df <- cbind(dcca_res2$distinct_score_1[,j], dcca_res2$common_score)
  df <- as.data.frame(df)
  colnames(df) <- c("y", paste0("x", 1:ncol(dcca_res2$common_score)))
  lm_res <- stats::lm(y ~ . , data = df)
  stats::residuals(lm_res)
})
alignment_1 <- 1 - sum((residual_mat_1)^2)/sum((dcca_res2$distinct_score_1)^2)

residual_mat_2 <- sapply(1:ncol(dcca_res2$distinct_score_2), function(j){
  df <- cbind(dcca_res2$distinct_score_2[,j], dcca_res2$common_score)
  df <- as.data.frame(df)
  colnames(df) <- c("y", paste0("x", 1:ncol(dcca_res2$common_score)))
  lm_res <- stats::lm(y ~ ., data = df)
  stats::residuals(lm_res)
})
alignment_2 <- 1 - sum((residual_mat_2)^2)/sum((dcca_res2$distinct_score_2)^2)

#######################################

n <- nrow(SNU@meta.data)
svd_1 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_1, 
                                       K = rank_1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_2, K = rank_2, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)
dimred_1 <- scale(dimred_1, center = T, scale = T)
tmp <- svd(dimred_1)
dimred_1 <- dimred_1/tmp$d[1]*sqrt(n)
l2_vec <- apply(dimred_1, 1, function(x){multiomicCCA:::.l2norm(x)})
dimred_1 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, dimred_1)
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*sqrt(n))
l2_vec <- apply(dimred_2, 1, function(x){multiomicCCA:::.l2norm(x)})
dimred_2 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, dimred_2)

set.seed(10)
common_umap <- Seurat::RunUMAP(cbind(dimred_1, dimred_2), 
                               metric = "euclidean",
                               reduction.key = "umapCommon_")
rownames(common_umap@cell.embeddings) <- rownames(SNU@meta.data)
SNU[["dcca_common"]] <- Seurat::CreateDimReducObject(common_umap@cell.embeddings, assay = "atac")

#######################################

svd_1 <- multiomicCCA:::.svd_truncated(dcca_decomp$distinct_mat_1, 
                                       K = rank_1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)
dimred_1 <- scale(dimred_1, center = T, scale = T)

set.seed(10)
distinct1_umap <- Seurat::RunUMAP(dimred_1, 
                                  metric = "cosine",
                                  reduction.key = "umapDistinct1_")
rownames(distinct1_umap@cell.embeddings) <- rownames(SNU@meta.data)
SNU[["dcca_distinct1"]] <- Seurat::CreateDimReducObject(distinct1_umap@cell.embeddings, assay = "atac")

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
rownames(distinct2_umap@cell.embeddings) <- rownames(SNU@meta.data)
SNU[["dcca_distinct2"]] <- Seurat::CreateDimReducObject(distinct2_umap@cell.embeddings, assay = "atac")

#######################################

anchor_name <- "umap.atac"
other_names <- c("umap.cna", "wnn.umap", "dcca_common", "dcca_distinct1", "dcca_distinct2")

for(umap_name in other_names){
  print(umap_name)
  u_mat1 <- SNU[[anchor_name]]@cell.embeddings
  u_mat2 <- SNU[[umap_name]]@cell.embeddings
  tmp <- svd(t(u_mat1) %*% u_mat2)
  rotation_mat <- tmp$u %*% t(tmp$v)
  tmp <- u_mat2 %*% t(rotation_mat)
  rownames(tmp) <- rownames(SNU@meta.data)
  colnames(tmp) <- colnames(SNU[[umap_name]]@cell.embeddings)
  SNU[[umap_name]]@cell.embeddings <- tmp
}

# plot according to clones
reduction_vec <- c(anchor_name, other_names)
group_vec <- c("clone")
main_vec <- c("(ATAC)", "(Copy number)", "(WNN)", 
              "(Tilted-CCA, Common)",  
              paste0("(T-CCA, Distinct 1, Alignment: ", round(alignment_1, 2),")"), 
              paste0("(T-CCA, Distinct 2, Alignment: ", round(alignment_2, 2),")"))
file_vec <- c("atac", "cna", "wnn", "tiltedcca-common", "tiltedcca-distinct1", "tiltedcca-distinct2")

for(i in 1:length(reduction_vec)){
  for(j in 1:length(group_vec)){
    plot1 <- Seurat::DimPlot(SNU, reduction = reduction_vec[i],
                             group.by = group_vec[j], label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU:\n", main_vec[i]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14f/Writeup14f_SNU_", file_vec[i], "_", group_vec[j], ".png"),
                    plot1, device = "png", width = 5, height = 5, units = "in")
  }
}

###################################################

reduction_vec <- c("dcca_common", "dcca_distinct1", "dcca_distinct2")
file_vec <- c("tiltedcca-common", "tiltedcca-distinct1", "tiltedcca-distinct2")

for(i in 1:length(reduction_vec)){
  p_list <- lapply(sort(unique(SNU@meta.data[,"clone"])), function(clone){
    p0 <- Seurat::DimPlot(SNU, reduction = reduction_vec[i],
                          cells.highlight = rownames(SNU@meta.data)[which(SNU[["clone"]] == clone)])
    p0 <- p0 + ggplot2::theme(legend.position="none") + ggplot2::ggtitle(paste0("Clone ", clone))
    p0
  })
  
  p <- cowplot::plot_grid(p_list[[1]], p_list[[2]], p_list[[3]], 
                          p_list[[4]], p_list[[5]], p_list[[6]])
  cowplot::save_plot(filename = paste0("../../../../out/figures/Writeup14f/Writeup14f_SNU_", file_vec[i], "_separate.png"), p, 
                     ncol = 3, nrow = 2, base_asp = 1, device = "png")
}

##################

png(file = "../../../../out/figures/Writeup14f/Writeup14f_SNU_cca.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot(NA,
     xlim = c(0,1), ylim = c(0, 1),
     xlab = "CCA value (across latent dimensions)",
     ylab = "Tilt percentage (0, tilt towards Modality 2)",
     main = paste0("SNU\nT-CCA summary (", ncol(dcca_res2$common_score), " out of ", 
                   max(ncol(dcca_res2$distinct_score_1), ncol(dcca_res2$distinct_score_2)), " dims.)"))
lines(c(-1e5,1e5), rep(0.5, 2), col = 2, lty = 2)
lines(c(-1e5,1e5), rep(dcca_res$tilt_perc[1], 2), col = 3, lwd = 2)
points(dcca_res2$cca_obj, dcca_res2$tilt_perc, pch = 16)
graphics.off()

###########

png("../../../../out/figures/Writeup14f/Writeup14f_SNU_dcca_scores.png",
    height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
multiomicCCA::plot_scores_heatmap.dcca(dcca_res2,
                                       membership_vec = as.factor(SNU$clone),
                                       log_scale = T, scaling_power = 4)
graphics.off()



