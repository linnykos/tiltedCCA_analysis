rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

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

tmp <- multiomicCCA:::form_snns(num_neigh = 60,
                                svd_1 = svd_1,
                                svd_2 = svd_2,
                                bool_intersect = T,
                                distance_func = "cosine",
                                min_deg = 30,
                                verbose = T)
snn_mat_1 <- tmp$snn_mat_1; snn_mat_2 <- tmp$snn_mat_2
basis_1 <- multiomicCCA:::compute_laplacian_basis(snn_mat_1, k = 50, verbose = T)
basis_2 <- multiomicCCA:::compute_laplacian_basis(snn_mat_2, k = 50, verbose = T)

#########

set.seed(10)
umap_1 <- Seurat::RunUMAP(basis_1,
                          metric = "euclidean")
rownames(umap_1@cell.embeddings) <- rownames(bm@meta.data)
set.seed(10)
umap_2 <- Seurat::RunUMAP(basis_2,
                          metric = "euclidean")
rownames(umap_2@cell.embeddings) <- rownames(bm@meta.data)

bm[["lapRnaUMAP"]] <- Seurat::CreateDimReducObject(umap_1@cell.embeddings,
                                                   assay = "RNA")
plot1 <- Seurat::DimPlot(bm, reduction = "lapRnaUMAP",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Cite-seq):\nLaplacian of SNN for RNA"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_bm25_rna_laplacian_snn_umap2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

bm[["lapAdtUMAP"]] <- Seurat::CreateDimReducObject(umap_2@cell.embeddings,
                                                   assay = "RNA")
plot1 <- Seurat::DimPlot(bm, reduction = "lapAdtUMAP",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Cite-seq):\nLaplacian of SNN for ADT"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_bm25_adt_laplacian_snn_umap2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

###################

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

# common neighbor
set.seed(10)
common_mat <- multiomicCCA:::common_neighborhood(snn_mat_1 = snn_mat_1, 
                                                 snn_mat_2 = snn_mat_2,
                                                 clustering_1 = clustering_1, 
                                                 clustering_2 = clustering_2,
                                                 num_neigh = 30,
                                                 verbose = T)
quantile(sparseMatrixStats::rowSums2(common_mat))
common_basis <- multiomicCCA:::compute_laplacian_basis(common_mat, k = 50, verbose = T)

set.seed(10)
common_umap <- Seurat::RunUMAP(common_basis,
                          metric = "euclidean")
rownames(common_umap@cell.embeddings) <- rownames(bm@meta.data)
bm[["lapComUMAP"]] <- Seurat::CreateDimReducObject(common_umap@cell.embeddings,
                                                   assay = "RNA")
plot1 <- Seurat::DimPlot(bm, reduction = "lapComUMAP",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Cite-seq):\nLaplacian of SNN for Common"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_bm25_common_laplacian_snn_umap2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


png("../../../../out/figures/Writeup14j/Writeup14j_citeseq_bm25_laplacian_scores.png",
    height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
multiomicCCA::plot_scores_heatmap.list(list(basis_1, basis_2, common_basis),
                                       main_vec = c("RNA basis", "ADT basis", "Common basis"),
                                       membership_vec = as.factor(bm$celltype.l2),
                                       log_scale = T, scaling_power = 4)
graphics.off()
tmp <- list(basis_1, basis_2, common_basis)
for(i in 1:length(tmp)){
  l2_vec <- apply(tmp[[i]], 2, multiomicCCA:::.l2norm)
  tmp[[i]] <- multiomicCCA:::.mult_mat_vec(tmp[[i]], 1/l2_vec)
}
png("../../../../out/figures/Writeup14j/Writeup14j_citeseq_bm25_laplacian_scores2.png",
    height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
multiomicCCA::plot_scores_heatmap.list(tmp,
                                       main_vec = c("RNA basis\n(Normalized)", "ADT basis\n(Normalized)", "Common basis\n(Normalized)"),
                                       membership_vec = as.factor(bm$celltype.l2),
                                       log_scale = T, scaling_power = 4)
graphics.off()

#############################

set.seed(10)
dcca_res <- multiomicCCA::dcca_factor(mat_1b, mat_2b, 
                                      dims_1 = 1:rank_1, dims_2 = 1:rank_2,
                                      center_1 = F, center_2 = F,
                                      scale_1 = F, scale_2 = F,
                                      discretization_gridsize = 11,
                                      fix_tilt_perc = F, 
                                      enforce_boundary = F,
                                      target_dimred = common_basis,
                                      verbose = T)
dcca_res2 <- multiomicCCA:::fine_tuning(dcca_res, verbose = T)
dcca_res2$tilt_perc

save.image("../../../../out/Writeup14j/Writeup14j_citeseq_bm25_dcca.RData")

