rm(list=ls())
load("../../../../out/Writeup14e/Writeup14e_citeseq_bm_dcca.RData")

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res)

#######

n <- nrow(bm@meta.data)
svd_1 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_1, 
                                       K = rank_1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = T, sd_vec = T,
                                       symmetric = F)
svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_2, K = rank_2, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = T, sd_vec = T,
                                       symmetric = F)
dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*sqrt(n))
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*sqrt(n))

set.seed(10)
common_umap <- Seurat::RunUMAP(cbind(dimred_1, dimred_2), 
                               metric = "cosine",
                               reduction.key = "umapCommon_")
rownames(common_umap@cell.embeddings) <- rownames(bm@meta.data)
bm[["dcca_common"]] <- Seurat::CreateDimReducObject(common_umap@cell.embeddings)

u_mat1 <- bm[["rna.umap"]]@cell.embeddings
u_mat2 <- bm[["dcca_common"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(bm@meta.data)
colnames(tmp) <- colnames(bm[["dcca_common"]]@cell.embeddings)
bm[["dcca_common"]]@cell.embeddings <- tmp


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
rownames(distinct1_umap@cell.embeddings) <- rownames(bm@meta.data)
bm[["dcca_distinct1"]] <- Seurat::CreateDimReducObject(distinct1_umap@cell.embeddings)

u_mat1 <- bm[["rna.umap"]]@cell.embeddings
u_mat2 <- bm[["dcca_distinct1"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(bm@meta.data)
colnames(tmp) <- colnames(bm[["dcca_distinct1"]]@cell.embeddings)
bm[["dcca_distinct1"]]@cell.embeddings <- tmp


##########

svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$distinct_mat_2, 
                                       K = rank_2, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = T, sd_vec = T,
                                       symmetric = F)
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d)

set.seed(10)
distinct2_umap <- Seurat::RunUMAP(dimred_2, 
                                  metric = "cosine",
                                  reduction.key = "umapDistinct2_")
rownames(distinct2_umap@cell.embeddings) <- rownames(bm@meta.data)
bm[["dcca_distinct2"]] <- Seurat::CreateDimReducObject(distinct2_umap@cell.embeddings)

u_mat1 <- bm[["rna.umap"]]@cell.embeddings
u_mat2 <- bm[["dcca_distinct2"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(bm@meta.data)
colnames(tmp) <- colnames(bm[["dcca_distinct2"]]@cell.embeddings)
bm[["dcca_distinct2"]]@cell.embeddings <- tmp

#################################################

u_mat1 <- bm[["rna.umap"]]@cell.embeddings
u_mat2 <- bm[["adt.umap"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(bm@meta.data)
colnames(tmp) <- colnames(bm[["adt.umap"]]@cell.embeddings)
bm[["adt.umap"]]@cell.embeddings <- tmp

###########################



plot1 <- Seurat::DimPlot(bm, reduction = "rna.umap", 
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow:\nCITE-Seq (RNA)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14e/Writeup14e_citeseq_bm_rna_umap.png",
                plot1, device = "png", width = 6, height = 5, units = "in")


plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap", 
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow:\nCITE-Seq (Protein)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14e/Writeup14e_citeseq_bm_atac_umap.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(bm, reduction = "dcca_common", 
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow:\nCITE-Seq (Tilted-CCA, Common)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14e/Writeup14e_citeseq_bm_tiltedcca_common_umap2.png",
                plot1, device = "png", width = 6, height = 5, units = "in")


plot1 <- Seurat::DimPlot(bm, reduction = "dcca_distinct1", 
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow:\nCITE-Seq (Tilted-CCA, Distinct 1)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14e/Writeup14e_citeseq_bm_tiltedcca_distinct1_umap2.png",
                plot1, device = "png", width = 6, height = 5, units = "in")


plot1 <- Seurat::DimPlot(bm, reduction = "dcca_distinct2", 
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow:\nCITE-Seq (Tilted-CCA, Distinct 2)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14e/Writeup14e_citeseq_bm_tiltedcca_distinct2_umap2.png",
                plot1, device = "png", width = 6, height = 5, units = "in")



