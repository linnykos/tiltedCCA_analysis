rm(list=ls())
load("../../../out/main/citeseq_pbmc224_preprocessed.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(1, ncol(pbmc))
keep_vec[which(pbmc$celltype.l2 == "Doublet")] <- 0
pbmc$keep <- keep_vec
pbmc <- subset(pbmc, keep == 1)

##########################################

Seurat::DefaultAssay(pbmc) <- "SCT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:40, reduction = "pca")
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)

Seurat::DefaultAssay(pbmc) <- "ADT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:25, reduction = "apca")
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)

plot1 <- Seurat::DimPlot(pbmc, reduction = "rna.umap",
                         group.by = "SCT_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster = F)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nSubset, RNA clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_rna-clustering.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")

plot1 <- Seurat::DimPlot(pbmc, reduction = "adt.umap",
                         group.by = "ADT_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster = F)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nSubset, ADT clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_adt-clustering.png"),
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
                                           num_metacells = 5000,
                                           verbose = 1)
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 30,
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

#############################################

svd_1 <- multiSVD_obj$svd_1
svd_2 <- multiSVD_obj$svd_2
consensus_pca <- tiltedCCA:::consensus_pca(mat_1 = NULL, mat_2 = NULL,
                                           dims_1 = NULL, dims_2 = NULL,
                                           dims_consensus = 1:max(ncol(svd_1$u), ncol(svd_2$u)),
                                           svd_1 = svd_1, svd_2 = svd_2,
                                           verbose = 1)

set.seed(10)
umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(pbmc)
pbmc[["consensusUMAP"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                        assay = "SCT",
                                                        key = "cUMAP")
source("pbmc_224antibody_colorPalette.R")

plot1 <- Seurat::DimPlot(pbmc, reduction = "consensusUMAP",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nConsensus PCA"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_consensusPCA-umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

