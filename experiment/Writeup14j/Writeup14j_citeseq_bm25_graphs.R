rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)
source("snn.R")

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
dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)
norm_vec <- apply(dimred_1, 1, multiomicCCA:::.l2norm)
dimred_1 <- multiomicCCA:::.mult_vec_mat(1/norm_vec, dimred_1)

svd_2 <- multiomicCCA:::.check_svd(svd_2, dims = dims_2)
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d)
norm_vec <- apply(dimred_2, 1, multiomicCCA:::.l2norm)
dimred_2 <- multiomicCCA:::.mult_vec_mat(1/norm_vec, dimred_2)

###########

snn_1 <- form_snn_mat(dimred_1, num_neigh = 60, 
                      bool_intersect = F,
                      min_deg = 30,
                      verbose = T)
basis_1 <- compute_laplacian_basis(snn_1, k = 50, verbose = T)
set.seed(10)
umap_1 <- Seurat::RunUMAP(basis_1,
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

snn_2 <- form_snn_mat(dimred_2, num_neigh = 60, 
                      bool_intersect = F,
                      min_deg = 30,
                      verbose = T)
basis_2 <- compute_laplacian_basis(snn_2, k = 50, verbose = T)
set.seed(10)
umap_2 <- Seurat::RunUMAP(basis_2,
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


####################################3

n <- ncol(bm)

Seurat::DefaultAssay(bm) <- "RNA"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:30)
bm <- Seurat::FindClusters(bm, resolution = 0.25)

Seurat::DefaultAssay(bm) <- "ADT"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:18, reduction = "apca")
bm <- Seurat::FindClusters(bm, resolution = 0.25)

clustering_1 <- factor(bm$RNA_snn_res.0.25)
clustering_2 <- factor(bm$ADT_snn_res.0.25)
if(any(table(clustering_1) == 0)) clustering_1 <- factor(clustering_1)
if(any(table(clustering_2) == 0)) clustering_2 <- factor(clustering_2)

snn_mat_1 <- snn_1
snn_mat_2 <- snn_2

snn_mat <- kl_selection(snn_mat_1 = snn_mat_1, 
                        snn_mat_2 = snn_mat_2,
                        clustering_1 = clustering_1, 
                        clustering_2 = clustering_2,
                        num_neigh = 30,
                        verbose = T)

common_basis <- compute_laplacian_basis(snn_mat, k = 50, verbose = T)
set.seed(10)
umap_common <- Seurat::RunUMAP(common_basis,
                          metric = "euclidean")
rownames(umap_common@cell.embeddings) <- rownames(bm@meta.data)
bm[["lapComUMAP"]] <- Seurat::CreateDimReducObject(umap_common@cell.embeddings,
                                                   assay = "RNA")
plot1 <- Seurat::DimPlot(bm, reduction = "lapComUMAP",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Cite-seq):\nLaplacian of SNN for Common"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_bm25_common_laplacian_snn_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

########################

