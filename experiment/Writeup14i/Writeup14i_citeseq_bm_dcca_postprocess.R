rm(list=ls())
library(Seurat); library(Signac)
load("../../../../out/Writeup14i/Writeup14i_citeseq_bm_dcca.RData")
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

n <- nrow(bm@meta.data)
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
rownames(common_umap@cell.embeddings) <- rownames(bm@meta.data)
bm[["dcca_common"]] <- Seurat::CreateDimReducObject(common_umap@cell.embeddings)

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
rownames(distinct1_umap@cell.embeddings) <- rownames(bm@meta.data)
bm[["dcca_distinct1"]] <- Seurat::CreateDimReducObject(distinct1_umap@cell.embeddings)

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
rownames(distinct2_umap@cell.embeddings) <- rownames(bm@meta.data)
bm[["dcca_distinct2"]] <- Seurat::CreateDimReducObject(distinct2_umap@cell.embeddings)


#######################################

anchor_name <- "rna.umap"
other_names <- c("adt.umap", "wnn.umap", "dcca_common", "dcca_distinct1", "dcca_distinct2")

for(umap_name in other_names){
  print(umap_name)
  u_mat1 <- bm[[anchor_name]]@cell.embeddings
  u_mat2 <- bm[[umap_name]]@cell.embeddings
  tmp <- svd(t(u_mat1) %*% u_mat2)
  rotation_mat <- tmp$u %*% t(tmp$v)
  tmp <- u_mat2 %*% t(rotation_mat)
  rownames(tmp) <- rownames(bm@meta.data)
  colnames(tmp) <- colnames(bm[[umap_name]]@cell.embeddings)
  bm[[umap_name]]@cell.embeddings <- tmp
}

# plot according to clones
reduction_vec <- c(anchor_name, other_names)
group_vec <- c("celltype.l2")
main_vec <- c("(RNA)", "(Protein)", "(WNN)", 
              "(Tilted-CCA, Common)", 
              paste0("(T-CCA, Distinct 1, Alignment: ", round(alignment_1, 2),")"), 
              paste0("(T-CCA, Distinct 2, Alignment: ", round(alignment_2, 2),")"))
file_vec <- c("rna", "adt", "wnn", "tiltedcca-common", "tiltedcca-distinct1", "tiltedcca-distinct2")

for(i in 1:length(reduction_vec)){
  for(j in 1:length(group_vec)){
    plot1 <- Seurat::DimPlot(bm, reduction = reduction_vec[i],
                             group.by = group_vec[j], label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Cite-seq):\n", main_vec[i]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14i/Writeup14i_citeseq_bm_", file_vec[i], ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}