rm(list=ls())
library(Seurat); library(Signac)
load("../../../../out/Writeup14f/Writeup14f_10x_mouseembryo_dcca.RData")
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
  lm_res <- stats::lm(y ~ . , data = df)
  stats::residuals(lm_res)
})
alignment_2 <- 1 - sum((residual_mat_2)^2)/sum((dcca_res2$distinct_score_2)^2)

#######################################

n <- nrow(mbrain2@meta.data)
svd_1 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_1, 
                                       K = rank_1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = T, sd_vec = T,
                                       symmetric = F)
svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_2, K = rank_2, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*sqrt(n))
l2_vec <- apply(dimred_1, 1, function(x){multiomicCCA:::.l2norm(x)})
dimred_1 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, dimred_1)
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u[,-1], svd_2$d[-1])
dimred_2 <- scale(dimred_2, center = T, scale = T)
tmp <- svd(dimred_2)
dimred_2 <- dimred_2/tmp$d[1]*sqrt(n)
l2_vec <- apply(dimred_2, 1, function(x){multiomicCCA:::.l2norm(x)})
dimred_2 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, dimred_2)

set.seed(10)
common_umap <- Seurat::RunUMAP(cbind(dimred_1, dimred_2), 
                               metric = "euclidean",
                               reduction.key = "umapCommon_")
rownames(common_umap@cell.embeddings) <- rownames(mbrain2@meta.data)
mbrain2[["dcca_common"]] <- Seurat::CreateDimReducObject(common_umap@cell.embeddings)

#######################################

svd_1 <- multiomicCCA:::.svd_truncated(dcca_decomp$distinct_mat_1, 
                                       K = rank_1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = T, sd_vec = T,
                                       symmetric = F)
dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)

set.seed(10)
distinct1_umap <- Seurat::RunUMAP(dimred_1, 
                                  metric = "cosine",
                                  reduction.key = "umapDistinct1_")
rownames(distinct1_umap@cell.embeddings) <- rownames(mbrain2@meta.data)
mbrain2[["dcca_distinct1"]] <- Seurat::CreateDimReducObject(distinct1_umap@cell.embeddings)

#######################################

svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$distinct_mat_2, 
                                       K = rank_2, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u[,-1], svd_2$d[-1])
dimred_2 <- scale(dimred_2, center = T, scale = T)

set.seed(10)
distinct2_umap <- Seurat::RunUMAP(dimred_2, 
                                  metric = "cosine",
                                  reduction.key = "umapDistinct2_")
rownames(distinct2_umap@cell.embeddings) <- rownames(mbrain2@meta.data)
mbrain2[["dcca_distinct2"]] <- Seurat::CreateDimReducObject(distinct2_umap@cell.embeddings)


#######################################

anchor_name <- "umap"
other_names <- c("umap.atac", "wnn.umap", "dcca_common", "dcca_distinct1", "dcca_distinct2")

for(umap_name in other_names){
  print(umap_name)
  u_mat1 <- mbrain2[[anchor_name]]@cell.embeddings
  u_mat2 <- mbrain2[[umap_name]]@cell.embeddings
  tmp <- svd(t(u_mat1) %*% u_mat2)
  rotation_mat <- tmp$u %*% t(tmp$v)
  tmp <- u_mat2 %*% t(rotation_mat)
  rownames(tmp) <- rownames(mbrain2@meta.data)
  colnames(tmp) <- colnames(mbrain2[[umap_name]]@cell.embeddings)
  mbrain2[[umap_name]]@cell.embeddings <- tmp
}

# plot according to clones
reduction_vec <- c(anchor_name, other_names)
group_vec <- c("label_Savercat")
main_vec <- c("(RNA)", "(ATAC)", "(WNN)", 
              "(Tilted-CCA, Common)",  
              paste0("(T-CCA, Distinct 1, Alignment: ", round(alignment_1, 2),")"), 
              paste0("(T-CCA, Distinct 2, Alignment: ", round(alignment_2, 2),")"))
file_vec <- c("rna", "atac", "wnn", "tiltedcca-common", "tiltedcca-distinct1", "tiltedcca-distinct2")

for(i in 1:length(reduction_vec)){
  for(j in 1:length(group_vec)){
    plot1 <- Seurat::DimPlot(mbrain2, reduction = reduction_vec[i],
                             group.by = group_vec[j], label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x):\n", main_vec[i]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14f/Writeup14f_10x_mouseembryo_", file_vec[i], "_", group_vec[j], ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}

########################################

png(file = "../../../../out/figures/Writeup14f/Writeup14f_10x_mouseembryo_cca.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot(NA,
     xlim = c(0,1), ylim = c(0, 1),
     xlab = "CCA value (across latent dimensions)",
     ylab = "Tilt percentage (0, tilt towards Modality 2)",
     main = paste0("Mouse Embryo (10x)\nT-CCA summary (", ncol(dcca_res2$common_score), " out of ", 
                   max(ncol(dcca_res2$distinct_score_1), ncol(dcca_res2$distinct_score_2)), " dims.)"))
lines(c(-1e5,1e5), rep(0.5, 2), col = 2, lty = 2)
lines(c(-1e5,1e5), rep(dcca_res$tilt_perc[1], 2), col = 3, lwd = 2)
points(dcca_res2$cca_obj, dcca_res2$tilt_perc, pch = 16)
graphics.off()


###########

png("../../../../out/figures/Writeup14f/Writeup14f_10x_mouseembryo_dcca_scores.png",
    height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
multiomicCCA::plot_scores_heatmap.dcca(dcca_res2,
                                       membership_vec = as.factor(mbrain2$label_Savercat),
                                       log_scale = T, scaling_power = 4)
graphics.off()