rm(list=ls())
load("../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_preprocessed.RData")
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
n <- ncol(pbmc)

membership_vec <- as.factor(pbmc@meta.data$celltype.l2)
n <- length(membership_vec)
sort(table(membership_vec))
set.seed(10)
idx <- multiomicCCA::construct_celltype_subsample(membership_vec, 
                                                  min_subsample_cell = 5000)
keep_vec <- rep(0, n)
names(keep_vec) <- rownames(pbmc@meta.data)
keep_vec[idx] <- 1
pbmc$keep <- keep_vec
pbmc <- subset(pbmc, keep == 1)
pbmc[["rna.umap"]] <- NULL
pbmc[["adt.umap"]] <- NULL
pbmc[["wnn.umap"]] <- NULL

##############

Seurat::DefaultAssay(pbmc) <- "SCT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:40)
set.seed(10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)
tab_vec <- table(pbmc$SCT_snn_res.0.25)
round(tab_vec/n, 3)
rm_idx <- names(tab_vec[tab_vec/n < 0.01])
pbmc$SCT_snn_res.0.25[pbmc$SCT_snn_res.0.25 %in% rm_idx] <- NA

Seurat::DefaultAssay(pbmc) <- "ADT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:50, reduction = "apca")
set.seed(10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)
tab_vec <- table(pbmc$ADT_snn_res.0.25)
round(tab_vec/n, 3)
rm_idx <- names(tab_vec[tab_vec/n < 0.01])
pbmc$ADT_snn_res.0.25[pbmc$ADT_snn_res.0.25 %in% rm_idx] <- NA

############

mat_1 <- t(pbmc[["SCT"]]@scale.data)
mat_2 <- t(pbmc[["ADT"]]@scale.data)

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

###################

dims_1 <- 1:40; dims_2 <- 1:25; nn <- 30
svd_1 <- multiomicCCA:::.svd_truncated(mat = mat_1, K = dims_1[2],
                                       symmetric = F, rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       K_full_rank = F)
svd_2 <- multiomicCCA:::.svd_truncated(mat = mat_2, K = dims_2[2],
                                       symmetric = F, rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       K_full_rank = F)
dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d)
metacell_clustering_1 = factor(pbmc$SCT_snn_res.0.25)
metacell_clustering_2 = factor(pbmc$ADT_snn_res.0.25)

min_res <- compute_min_subspace(dimred_1, dimred_2,
                                metacell_clustering_1 = metacell_clustering_1,
                                metacell_clustering_2 = metacell_clustering_2,
                                binarize = F,
                                num_neigh = nn,
                                verbose = T)

##################################

svd_tmp <- multiomicCCA:::.svd_truncated(mat = min_res$min_mat, 
                                         K = 25,
                                         symmetric = F, rescale = F,
                                         mean_vec = T, sd_vec = F,
                                         K_full_rank = F)
set.seed(10)
umap_res <- Seurat::RunUMAP(multiomicCCA:::.mult_mat_vec(svd_tmp$u, svd_tmp$d), 
                            metric = "euclidean",
                            assay = "SCT",
                            reduction.key = "umap1_")
graph_obj <- SeuratObject::as.Graph(min_res$min_mat)
set.seed(10)
umap_res2 <- Seurat::RunUMAP(graph_obj, 
                            metric = "euclidean",
                            assay = "SCT",
                            reduction.key = "umap2_")
pbmc[["umap1"]] <- Seurat::CreateDimReducObject(umap_res@cell.embeddings,
                                                   assay = "SCT")
pbmc[["umap2"]] <- Seurat::CreateDimReducObject(umap_res2@cell.embeddings,
                                                assay = "SCT")

plot1 <- Seurat::DimPlot(pbmc, reduction = "umap1",
                         group.by = "celltype.l2",
                         label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq, 224):\nMin-dist (UMAP-PCA)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot2 <- Seurat::DimPlot(pbmc, reduction = "umap2",
                         group.by = "celltype.l2",
                         label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq, 224):\nMin-dist (UMAP-graph)"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
p <- cowplot::plot_grid(plot1, plot2)
cowplot::save_plot(filename = paste0("../../../../out/figures/Writeup14h/Writeup14h_citeseq_pbmc_mindist_non-binarize.png"), p, 
                   ncol = 2, nrow = 1, base_asp = 1.2,
                   base_height = 6, device = "png")

###############3

# binary version
tmp_mat <- min_res$min_mat; tmp_mat@x <- rep(1, length(tmp_mat@x))
svd_tmp <- multiomicCCA:::.svd_truncated(mat = tmp_mat, 
                                         K = 25,
                                         symmetric = F, rescale = F,
                                         mean_vec = T, sd_vec = F,
                                         K_full_rank = F)
set.seed(10)
umap_res <- Seurat::RunUMAP(multiomicCCA:::.mult_mat_vec(svd_tmp$u, svd_tmp$d), 
                            metric = "euclidean",
                            assay = "SCT",
                            reduction.key = "umap1_")
graph_obj <- SeuratObject::as.Graph(tmp_mat)
set.seed(10)
umap_res2 <- Seurat::RunUMAP(graph_obj, 
                             metric = "euclidean",
                             assay = "SCT",
                             reduction.key = "umap2_")
pbmc[["umap1"]] <- Seurat::CreateDimReducObject(umap_res@cell.embeddings,
                                                assay = "SCT")
pbmc[["umap2"]] <- Seurat::CreateDimReducObject(umap_res2@cell.embeddings,
                                                assay = "SCT")

plot1 <- Seurat::DimPlot(pbmc, reduction = "umap1",
                         group.by = "celltype.l2",
                         label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq, 224):\nMin-dist (UMAP-PCA, bin)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot2 <- Seurat::DimPlot(pbmc, reduction = "umap2",
                         group.by = "celltype.l2",
                         label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq, 224):\nMin-dist (UMAP-graph, bin)"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
p <- cowplot::plot_grid(plot1, plot2)
cowplot::save_plot(filename = paste0("../../../../out/figures/Writeup14h/Writeup14h_citeseq_pbmc_mindist_binarize.png"), p, 
                   ncol = 2, nrow = 1, base_asp = 1.2,
                   base_height = 6, device = "png")



# dims_1 <- 1:40; dims_2 <- 1:25; nn <- 30
# set.seed(10)
# 
# max_mat <- zongming_embedding(mat_1, mat_2,
#                               dims_1, dims_2, 
#                               nn)
# 
# ########################3
# 
# 
# graph_obj <- SeuratObject::as.Graph(max_mat)
# 
# set.seed(10)
# umap_res <- Seurat::RunUMAP(graph_obj, 
#                             metric = "euclidean",
#                             assay = "SCT",
#                             reduction.key = "umapMax_")
# rownames(umap_res@cell.embeddings) <- rownames(mat_1)
# pbmc[["max.umap"]] <- Seurat::CreateDimReducObject(umap_res@cell.embeddings,
#                                                    assay = "SCT")
# 
# plot1 <- Seurat::DimPlot(pbmc, reduction = "max.umap", 
#                          group.by = "celltype.l2",
#                          label = TRUE, 
#                          repel = TRUE, label.size = 2.5)
# plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq, 224):\nMaximum distance"))
# plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
# ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14h/Writeup14h_citeseq_pbmc_v3_maxdist.png"),
#                 plot1, device = "png", width = 6, height = 5, units = "in")




