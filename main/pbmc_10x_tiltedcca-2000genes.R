rm(list=ls())
load("../../../out/main/10x_pbmc_preprocessed.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

Seurat::DefaultAssay(pbmc) <- "SCT"
pbmc[["SCT"]]@var.features <- pbmc[["SCT"]]@var.features[1:2000]
pbmc[["SCT"]]@scale.data <- pbmc[["SCT"]]@scale.data[pbmc[["SCT"]]@var.features,]

pbmc[["pca"]] <- NULL
pbmc[["umap.rna"]] <- NULL
set.seed(10)
pbmc <- Seurat::RunPCA(pbmc, verbose = FALSE)
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:50, assay = 'SCT',
                        reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

###################

Seurat::DefaultAssay(pbmc) <- "SCT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:50, reduction = "pca")
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)

Seurat::DefaultAssay(pbmc) <- "ATAC"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 2:50, reduction = "lsi")
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)

plot1 <- Seurat::DimPlot(pbmc, reduction = "umap.rna",
                         group.by = "SCT_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (10x, RNA+ATAC)\nRNA clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_pbmc-2000genes_rna-clustering.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")

plot1 <- Seurat::DimPlot(pbmc, reduction = "umap.atac",
                         group.by = "ATAC_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (10x, RNA+ATAC)\nATAC clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_pbmc-2000genes_atac-clustering.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")

#########

Seurat::DefaultAssay(pbmc) <- "SCT"
mat_1 <- Matrix::t(pbmc[["SCT"]]@data[Seurat::VariableFeatures(object = pbmc),])
Seurat::DefaultAssay(pbmc) <- "ATAC"
mat_2 <- Matrix::t(pbmc[["ATAC"]]@data[Seurat::VariableFeatures(object = pbmc),])

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
                                            dims_1 = 1:50, dims_2 = 2:50,
                                            center_1 = T, center_2 = F,
                                            normalize_row = T,
                                            normalize_singular_value = T,
                                            recenter_1 = F, recenter_2 = T,
                                            rescale_1 = F, rescale_2 = T,
                                            scale_1 = T, scale_2 = F,
                                            verbose = 1)
multiSVD_obj <- tiltedCCA:::form_metacells(input_obj = multiSVD_obj,
                                           large_clustering_1 = as.factor(pbmc$SCT_snn_res.0.25), 
                                           large_clustering_2 = as.factor(pbmc$ATAC_snn_res.0.25), 
                                           num_metacells = 5000,
                                           verbose = 1)
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 50,
                                         num_neigh = 60,
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
                                                    bool_modality_2_full = F)

save(multiSVD_obj, pbmc,
     date_of_run, session_info,
     file = "../../../out/main/10x_pbmc_tiltedcca-2000genes.RData")

#################################

set.seed(10)
pbmc[["common_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                      what = "common",
                                                      aligned_umap_assay = "umap.rna",
                                                      seurat_obj = pbmc,
                                                      seurat_assay = "SCT",
                                                      verbose = 1)
set.seed(10)
pbmc[["distinct1_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                         what = "distinct_1",
                                                         aligned_umap_assay = "umap.rna",
                                                         seurat_obj = pbmc,
                                                         seurat_assay = "SCT",
                                                         verbose = 1)
set.seed(10)
pbmc[["distinct2_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                         what = "distinct_2",
                                                         aligned_umap_assay = "umap.rna",
                                                         seurat_obj = pbmc,
                                                         seurat_assay = "SCT",
                                                         verbose = 1)

save(multiSVD_obj, pbmc,
     date_of_run, session_info,
     file = "../../../out/main/10x_pbmc_tiltedcca-2000genes.RData")
