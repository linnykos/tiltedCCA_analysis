rm(list=ls())
load("../../../out/main/citeseq_pbmc224_preprocessed.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

membership_vec <- as.factor(pbmc$celltype.l2)
n <- length(membership_vec)
set.seed(10)
idx <- tiltedCCA::construct_celltype_subsample(membership_vec, 
                                               min_subsample_cell = 5000)
keep_vec <- rep(0, n)
names(keep_vec) <- rownames(pbmc@meta.data)
keep_vec[idx] <- 1
pbmc$keep <- keep_vec
pbmc <- subset(pbmc, keep == 1)
pbmc[["rna.umap"]] <- NULL
pbmc[["adt.umap"]] <- NULL
pbmc[["wnn.umap"]] <- NULL

##########################################

Seurat::DefaultAssay(pbmc) <- "SCT"
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'pca', dims = 1:40, assay = 'SCT',
                        reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:40, reduction = "pca")
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)

Seurat::DefaultAssay(pbmc) <- "ADT"
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'apca', dims = 1:25, assay = 'ADT',
                        reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:25, reduction = "apca")
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)


plot1 <- Seurat::DimPlot(pbmc, reduction = "rna.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nSubset, RNA UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_subset_rna-umap.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")
plot1 <- Seurat::DimPlot(pbmc, reduction = "rna.umap",
                         group.by = "SCT_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nSubset, RNA clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_subset_rna-clustering.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")

plot1 <- Seurat::DimPlot(pbmc, reduction = "adt.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nSubset, ADT UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_subset_adt-umap.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")
plot1 <- Seurat::DimPlot(pbmc, reduction = "adt.umap",
                         group.by = "ADT_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nSubset, ADT clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_subset_adt-clustering.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")

#########

Seurat::DefaultAssay(pbmc) <- "SCT"
mat_1 <- Matrix::t(pbmc[["SCT"]]@scale.data)
Seurat::DefaultAssay(pbmc) <- "ADT"
mat_2 <- Matrix::t(pbmc[["ADT"]]@scale.data)

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

#################################

set.seed(10)
multiSVD_obj <- tiltedCCA:::create_multiSVD(mat_1 = mat_1b, mat_2 = mat_2b,
                                            dims_1 = 1:40, dims_2 = 1:25,
                                            center_1 = T, center_2 = T,
                                            normalize_row = T,
                                            normalize_singular_value = T,
                                            recenter_1 = F, recenter_2 = F,
                                            rescale_1 = F, rescale_2 = F,
                                            scale_1 = T, scale_2 = T,
                                            verbose = 1)
multiSVD_obj <- tiltedCCA:::form_metacells(input_obj = multiSVD_obj,
                                           large_clustering_1 = as.factor(pbmc$SCT_snn_res.0.25), 
                                           large_clustering_2 = as.factor(pbmc$ADT_snn_res.0.25), 
                                           num_metacells = NULL,
                                           verbose = 1)
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 50,
                                         num_neigh = 60,
                                         bool_cosine = T,
                                         bool_intersect = F,
                                         min_deg = 30,
                                         verbose = 2)
tmp <- pbmc; tmp_mat <- multiSVD_obj$laplacian_list$common_laplacian
colnames(tmp_mat) <- paste0("tmp_", 1:ncol(tmp_mat))
set.seed(10); tmp_umap <- Seurat::RunUMAP(tmp_mat)
tmp[["common_laplacian"]] <- Seurat::CreateDimReducObject(tmp_umap@cell.embeddings,
                                                          assay = "SCT")
plot1 <- Seurat::DimPlot(tmp, reduction = "common_laplacian",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nCommon Laplacian"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_umap_common-laplacian.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

multiSVD_obj <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj,
                                      verbose = 1)
multiSVD_obj <- tiltedCCA:::fine_tuning(input_obj = multiSVD_obj,
                                        verbose = 1)
multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = T,
                                                    bool_modality_2_full = T)

save(multiSVD_obj, pbmc,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_pbmc224_tiltedcca.RData")

#################################

set.seed(10)
pbmc[["common_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                      what = "common",
                                                      aligned_umap_assay = "rna.umap",
                                                      seurat_obj = pbmc,
                                                      seurat_assay = "SCT",
                                                      verbose = 1)
set.seed(10)
pbmc[["distinct1_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                         what = "distinct_1",
                                                         aligned_umap_assay = "rna.umap",
                                                         seurat_obj = pbmc,
                                                         seurat_assay = "SCT",
                                                         verbose = 1)
set.seed(10)
pbmc[["distinct2_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                         what = "distinct_2",
                                                         aligned_umap_assay = "rna.umap",
                                                         seurat_obj = pbmc,
                                                         seurat_assay = "SCT",
                                                         verbose = 1)

save(multiSVD_obj, pbmc,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_pbmc224_tiltedcca.RData")
