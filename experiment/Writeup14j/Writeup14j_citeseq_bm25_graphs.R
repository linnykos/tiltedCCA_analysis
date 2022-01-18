rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)
library(multiomicCCA)

rank_1 <- 30; rank_2 <- 18; num_neigh <- 30
dims_1 <- 1:rank_1; dims_2 <- 1:rank_2

mat_1 <- t(bm[["RNA"]]@scale.data)
mat_2 <- t(bm[["ADT"]]@scale.data)
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

svd_1 <- multiomicCCA:::.svd_truncated(mat_1, K = rank_1, symmetric = F, rescale = F, 
                        mean_vec = F, sd_vec = F, K_full_rank = F)
svd_2 <- multiomicCCA:::.svd_truncated(mat_2, K = rank_2, symmetric = F, rescale = F, 
                        mean_vec = F, sd_vec = F, K_full_rank = F)

svd_1 <- multiomicCCA:::.check_svd(svd_1, dims = dims_1)
svd_2 <- multiomicCCA:::.check_svd(svd_2, dims = dims_2)

###########

tmp_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)
norm_vec <- apply(tmp_1, 1, multiomicCCA:::.l2norm)
tmp_1 <- multiomicCCA:::.mult_vec_mat(1/norm_vec, tmp_1)
snn_mat_1 <- multiomicCCA:::.form_snn_mat(bool_intersect = T,
                           mat = tmp_1, 
                           num_neigh = 30)
tmp_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d)
norm_vec <- apply(tmp_2, 1, multiomicCCA:::.l2norm)
tmp_2 <- multiomicCCA:::.mult_vec_mat(1/norm_vec, tmp_2)
snn_mat_2 <- multiomicCCA:::.form_snn_mat(bool_intersect = T,
                                          mat = tmp_2, 
                                          num_neigh = 30) ## [[Maybe this should be based on cosine distance??]]

deg_vec_1 <- sparseMatrixStats::rowSums2(snn_mat_1)
quantile(deg_vec_1)
length(which(deg_vec_1 == 0))
deg_vec_1[deg_vec_1 == 0] <- 1
diag_mat <- Matrix::Diagonal(x = 1/sqrt(deg_vec_1))
lap_mat_1 <-  diag_mat %*% snn_mat_1 %*% diag_mat

deg_vec_1 <- sparseMatrixStats::rowSums2(lap_mat_1)
deg_vec_1[deg_vec_1 == 0] <- 1
diag_mat <- Matrix::Diagonal(x = 1/deg_vec_1)
lap_mat_1 <-  diag_mat %*% lap_mat_1
eigen_1 <- irlba::partial_eigen(lap_mat_1, n = 50, symmetric = F)
tmp <- multiomicCCA:::.mult_mat_vec(eigen_1$vectors, eigen_1$values)
set.seed(10)
umap_1 <- Seurat::RunUMAP(tmp,
                          metric = "euclidean")
rownames(umap_1@cell.embeddings) <- rownames(bm@meta.data)
bm[["lapRnaUMAP"]] <- Seurat::CreateDimReducObject(umap_1@cell.embeddings,
                                                   assay = "RNA")
plot1 <- Seurat::DimPlot(bm, reduction = "lapRnaUMAP",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Cite-seq):\nLaplacian of SNN for RNA"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_bm25_rna_laplacian_snn_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

###############################

deg_vec_2 <- sparseMatrixStats::rowSums2(snn_mat_2)
quantile(deg_vec_2)
length(which(deg_vec_2 == 0))
deg_vec_2[deg_vec_2 == 0] <- 1
diag_mat <- Matrix::Diagonal(x = 1/sqrt(deg_vec_2))
lap_mat_2 <-  diag_mat %*% snn_mat_2 %*% diag_mat

deg_vec_2 <- sparseMatrixStats::rowSums2(lap_mat_2)
deg_vec_2[deg_vec_2 == 0] <- 1
diag_mat <- Matrix::Diagonal(x = 1/deg_vec_2)
lap_mat_2 <-  diag_mat %*% lap_mat_2
eigen_2 <- irlba::partial_eigen(lap_mat_2, n = 50, symmetric = F)
tmp <- multiomicCCA:::.mult_mat_vec(eigen_2$vectors, eigen_2$values)
set.seed(10)
umap_2 <- Seurat::RunUMAP(tmp,
                          metric = "euclidean")
rownames(umap_2@cell.embeddings) <- rownames(bm@meta.data)
bm[["lapAdtUMAP"]] <- Seurat::CreateDimReducObject(umap_2@cell.embeddings,
                                                   assay = "ADT")
plot1 <- Seurat::DimPlot(bm, reduction = "lapAdtUMAP",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Cite-seq):\nLaplacian of SNN for ADT"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_bm25_adt_laplacian_snn_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")





