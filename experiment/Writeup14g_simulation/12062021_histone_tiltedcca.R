rm(list=ls())
load("../../out/simulation/Writeup14g_data.RData")

mat_1 <- seurat_obj2[["pca"]]@cell.embeddings
mat_2 <- seurat_obj2[["lsi"]]@cell.embeddings
rank_1 <- 20; rank_2 <- 20; nn <- 30
set.seed(10)
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, 
                                      dims_1 = 1:rank_1, dims_2 = 1:rank_2,
                                      center_1 = F, center_2 = F,
                                      scale_1 = F, scale_2 = F,
                                      num_neigh = nn, 
                                      metacell_clustering_1 = NA,
                                      metacell_clustering_2 = NA,
                                      fix_tilt_perc = F, verbose = T)
dcca_res2 <- multiomicCCA:::fine_tuning(dcca_res, verbose = T)
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2, rank_c = 20)

set.seed(10)
dcca_common_umap <- Seurat::RunUMAP(cbind(dcca_decomp$common_mat_1, dcca_decomp$common_mat_2), 
                                    metric = "euclidean",
                                    reduction.key = "dccaCommon_")
seurat_obj2[["dcca_common"]] <- Seurat::CreateDimReducObject(dcca_common_umap@cell.embeddings)

u_mat1 <- seurat_obj2[["umap1"]]@cell.embeddings
u_mat2 <- seurat_obj2[["dcca_common"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("dccaCommon_1", "dccaCommon_2")
seurat_obj2[["dcca_common"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "dcca_common", 
                         cols = color_vec, 
                         group.by = "celltype_custom",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Tilted-CCA\n(Common)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/tiltedcca_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


set.seed(10)
dcca_distinct1_umap <- Seurat::RunUMAP(dcca_decomp$distinct_mat_1, 
                                       metric = "euclidean",
                                       reduction.key = "dccaDistinct1_")
seurat_obj2[["dcca_distinct1"]] <- Seurat::CreateDimReducObject(dcca_distinct1_umap@cell.embeddings)

u_mat1 <- seurat_obj2[["umap1"]]@cell.embeddings
u_mat2 <- seurat_obj2[["dcca_distinct1"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("dccaDistinct1_1", "dccaDistinct1_2")
seurat_obj2[["dcca_distinct1"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "dcca_distinct1", 
                         cols = color_vec, 
                         group.by = "celltype_custom",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Tilted-CCA\n(Distinct 1)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/tiltedcca_distinct1.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


set.seed(10)
dcca_distinct2_umap <- Seurat::RunUMAP(dcca_decomp$distinct_mat_2, 
                                       metric = "euclidean",
                                       reduction.key = "dccaDistinct2_")
seurat_obj2[["dcca_distinct2"]] <- Seurat::CreateDimReducObject(dcca_distinct2_umap@cell.embeddings)

u_mat1 <- seurat_obj2[["umap1"]]@cell.embeddings
u_mat2 <- seurat_obj2[["dcca_distinct2"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("dccaDistinct2_1", "dccaDistinct2_2")
seurat_obj2[["dcca_distinct2"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "dcca_distinct2", 
                         cols = color_vec, 
                         group.by = "celltype_custom",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Tilted-CCA\n(Distinct 2)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/tiltedcca_distinct2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")
