rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(multiomicCCA)

load("../../../../out/Writeup14d/Writeup14d_10x_mouseembryo_seurat.RData")
date_of_run <- Sys.time(); session_info <- devtools::session_info()

Seurat::DefaultAssay(mbrain2) <- "SCT"
mat_1 <- Matrix::t(mbrain2[["SCT"]]@data[Seurat::VariableFeatures(object = mbrain2),])
Seurat::DefaultAssay(mbrain2) <- "ATAC"
mat_2 <- Matrix::t(mbrain2[["ATAC"]]@data[Seurat::VariableFeatures(object = mbrain2),])

mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}

mat_2b <- mat_2
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}

##################

rank_1 <- 30; rank_2 <- 50; nn <- 15
svd_1 <- multiomicCCA:::.svd_truncated(mat_1b, K = rank_1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = T, sd_vec = T,
                                       symmetric = F)
tmp_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)
snn_mat_1 <- multiomicCCA:::.form_snn_mat(bool_intersect = T,
                           mat = tmp_1, 
                           num_neigh = nn)
metacell_clustering_1 <- lapply(1:nrow(snn_mat_1), function(i){
  which(snn_mat_1[i,] != 0)
})

svd_2 <- multiomicCCA:::.svd_truncated(mat_2b, K = rank_2, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
tmp_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u[,-1], svd_2$d[-1])
snn_mat_2 <- multiomicCCA:::.form_snn_mat(bool_intersect = T,
                           mat = tmp_2, 
                           num_neigh = nn)
metacell_clustering_2 <- lapply(1:nrow(snn_mat_2), function(i){
  which(snn_mat_2[i,] != 0)
})



##################

set.seed(10)
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, 
                                      dims_1 = 1:rank_1, dims_2 = 2:rank_2,
                                      center_1 = T, center_2 = F,
                                      scale_1 = T, scale_2 = F,
                                      num_neigh = nn, 
                                      metacell_clustering_1 = metacell_clustering_1,
                                      metacell_clustering_2 = metacell_clustering_2,
                                      fix_tilt_perc = F, verbose = T)
dcca_res$cca_obj
dcca_res$df_percentage
dcca_res$tilt_perc

dcca_res2 <- multiomicCCA:::fine_tuning(dcca_res, verbose = T)
dcca_res2$tilt_perc

save(mbrain2, dcca_res, dcca_res2, metacell_clustering_1, metacell_clustering_2,
     rank_1, rank_2, nn, date_of_run, session_info,
     file = "../../../../out/Writeup14e/Writeup14e_10x_mouseembryo_dcca.RData")

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2)

#####

n <- nrow(mbrain2@meta.data)
svd_1 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_1, 
                                       K = rank_1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = T, sd_vec = T,
                                       symmetric = F)
svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_2, K = rank_2-1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*sqrt(n))
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d)
dimred_2 <- scale(dimred_2, center = T, scale = F)
dimred_2 <- dimred_2/svd(dimred_2)$d[1]*sqrt(n)

set.seed(10)
common_umap <- Seurat::RunUMAP(cbind(dimred_1, dimred_2), 
                               metric = "cosine",
                               reduction.key = "umapCommon_")
rownames(common_umap@cell.embeddings) <- rownames(mbrain2@meta.data)
mbrain2[["dcca_common"]] <- Seurat::CreateDimReducObject(common_umap@cell.embeddings)

u_mat1 <- mbrain2[["umap"]]@cell.embeddings
u_mat2 <- mbrain2[["dcca_common"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(mbrain2@meta.data)
colnames(tmp) <- colnames(mbrain2[["dcca_common"]]@cell.embeddings)
mbrain2[["dcca_common"]]@cell.embeddings <- tmp

##########

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

u_mat1 <- mbrain2[["umap"]]@cell.embeddings
u_mat2 <- mbrain2[["dcca_distinct1"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(mbrain2@meta.data)
colnames(tmp) <- colnames(mbrain2[["dcca_distinct1"]]@cell.embeddings)
mbrain2[["dcca_distinct1"]]@cell.embeddings <- tmp


##########

svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$distinct_mat_2, 
                                       K = rank_2-1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = T, sd_vec = T,
                                       symmetric = F)
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d)

set.seed(10)
distinct2_umap <- Seurat::RunUMAP(dimred_2, 
                                  metric = "cosine",
                                  reduction.key = "umapDistinct2_")
rownames(distinct2_umap@cell.embeddings) <- rownames(mbrain2@meta.data)
mbrain2[["dcca_distinct2"]] <- Seurat::CreateDimReducObject(distinct2_umap@cell.embeddings)

u_mat1 <- mbrain2[["umap"]]@cell.embeddings
u_mat2 <- mbrain2[["dcca_distinct2"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(mbrain2@meta.data)
colnames(tmp) <- colnames(mbrain2[["dcca_distinct2"]]@cell.embeddings)
mbrain2[["dcca_distinct2"]]@cell.embeddings <- tmp



###############################################

dimred_atac <- multiomicCCA:::.mult_mat_vec(dcca_res$svd_2$u, dcca_res$svd_2$d)
dimred_atac <- scale(dimred_atac, center = T, scale = F)

set.seed(10)
atac_umap <- Seurat::RunUMAP(dimred_atac, 
                             metric = "cosine",
                             reduction.key = "umapATAC_")

png(file = "../../../../out/figures/Writeup14e/test.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot(atac_umap@cell.embeddings[,1], atac_umap@cell.embeddings[,2], asp = T,
     col = as.factor(mbrain2$label_Savercat))
graphics.off()

#################################################

u_mat1 <- mbrain2[["umap"]]@cell.embeddings
u_mat2 <- mbrain2[["umap.atac"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(mbrain2@meta.data)
colnames(tmp) <- colnames(mbrain2[["umap.atac"]]@cell.embeddings)
mbrain2[["umap.atac"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(mbrain2, reduction = "umap", 
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo: 10x (RNA)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14e/Writeup14e_10x_mouseembryo_rna_umap.png",
                plot1, device = "png", width = 6.5, height = 5, units = "in")


plot1 <- Seurat::DimPlot(mbrain2, reduction = "umap.atac", 
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo: 10x (ATAC)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14e/Writeup14e_10x_mouseembryo_atac_umap.png",
                plot1, device = "png", width = 6.5, height = 5, units = "in")

plot1 <- Seurat::DimPlot(mbrain2, reduction = "dcca_common", 
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo: 10x\n(Tilted-CCA, Common)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14e/Writeup14e_10x_mouseembryo_tiltedcca_common_umap.png",
                plot1, device = "png", width = 6.5, height = 5, units = "in")


plot1 <- Seurat::DimPlot(mbrain2, reduction = "dcca_distinct1", 
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo: 10x\n(Tilted-CCA, Distinct 1)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14e/Writeup14e_10x_mouseembryo_tiltedcca_distinct1_umap.png",
                plot1, device = "png", width = 6.5, height = 5, units = "in")


plot1 <- Seurat::DimPlot(mbrain2, reduction = "dcca_distinct2", 
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo: 10x\n(Tilted-CCA, Distinct 2)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14e/Writeup14e_10x_mouseembryo_tiltedcca_distinct2_umap.png",
                plot1, device = "png", width = 6.5, height = 5, units = "in")


