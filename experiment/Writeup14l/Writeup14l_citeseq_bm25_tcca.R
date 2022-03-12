rm(list=ls())
load("../../../../out/Writeup14k/Writeup14k_citeseq_bm25_preprocessed.RData")

library(Seurat)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
n <- ncol(bm)

Seurat::DefaultAssay(bm) <- "RNA"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:30)
bm <- Seurat::FindClusters(bm, resolution = 0.25)

Seurat::DefaultAssay(bm) <- "ADT"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:18, reduction = "apca")
bm <- Seurat::FindClusters(bm, resolution = 0.25)

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

#########################

set.seed(10)
multiSVD_obj <- tiltedCCA:::create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                            dims_1 = 1:30, dims_2 = 1:18,
                                            center_1 = T, center_2 = T,
                                            normalize_row = T,
                                            normalize_singular_value = T,
                                            recenter_1 = F, recenter_2 = F,
                                            rescale_1 = F, rescale_2 = F,
                                            scale_1 = T, scale_2 = T,
                                            verbose = 1)
multiSVD_obj <- tiltedCCA:::form_metacells(input_obj = multiSVD_obj,
                                           large_clustering_1 = as.factor(bm$RNA_snn_res.0.25), 
                                           large_clustering_2 = as.factor(bm$ADT_snn_res.0.25), 
                                           num_metacells = 5000,
                                           verbose = 1)
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 20,
                                         num_neigh = 30,
                                         bool_cosine = T,
                                         bool_intersect = F,
                                         min_deg = 30,
                                         verbose = 2)
multiSVD_obj <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj,
                                      verbose = 1)
multiSVD_obj <- tiltedCCA:::fine_tuning(input_obj = multiSVD_obj,
                                        verbose = 1)
multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1)

#################################

set.seed(10)
bm[["common_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                    what = "common",
                                                    aligned_umap_assay = "rna.umap",
                                                    seurat_obj = bm,
                                                    seurat_assay = "RNA",
                                                    verbose = 1)
set.seed(10)
bm[["distinct1_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                       what = "distinct_1",
                                                       aligned_umap_assay = "rna.umap",
                                                       seurat_obj = bm,
                                                       seurat_assay = "RNA",
                                                       verbose = 1)
set.seed(10)
bm[["distinct2_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                       what = "distinct_2",
                                                       aligned_umap_assay = "rna.umap",
                                                       seurat_obj = bm,
                                                       seurat_assay = "RNA",
                                                       verbose = 1)


plot1 <-Seurat::DimPlot(bm, reduction = "common_tcca",
                        group.by = "celltype.l2", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (CITE-Seq, RNA+ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_citeseq_bm25_umap_common-tcca.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "distinct1_tcca",
                        group.by = "celltype.l2", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (CITE-Seq, RNA+ADT)\nDistinct 1 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_citeseq_bm25_umap_distinct1-tcca.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "distinct2_tcca",
                        group.by = "celltype.l2", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (CITE-Seq, RNA+ADT)\nDistinct 2 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_citeseq_bm25_umap_distinct2-tcca.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")


