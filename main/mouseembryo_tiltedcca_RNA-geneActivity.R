rm(list=ls())
load("../../../out/main/10x_mouseembryo_preprocessed.RData")

library(Seurat); library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

Seurat::DefaultAssay(mbrain) <- "SCT"
mat_1 <- Matrix::t(mbrain[["SCT"]]@data[Seurat::VariableFeatures(object = mbrain),])
Seurat::DefaultAssay(mbrain) <- "geneActivity"
mat_2 <- Matrix::t(mbrain[["geneActivity"]]@data[Seurat::VariableFeatures(object = mbrain),])

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
                                            center_1 = T, center_2 = T,
                                            normalize_row = T,
                                            normalize_singular_value = T,
                                            recenter_1 = F, recenter_2 = F,
                                            rescale_1 = F, rescale_2 = F,
                                            scale_1 = T, scale_2 = T,
                                            verbose = 1)
multiSVD_obj <- tiltedCCA:::form_metacells(input_obj = multiSVD_obj,
                                           large_clustering_1 = NULL, 
                                           large_clustering_2 = NULL, 
                                           num_metacells = NULL,
                                           verbose = 1)
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 20,
                                         num_neigh = 30,
                                         bool_cosine = T,
                                         bool_intersect = F,
                                         min_deg = 30,
                                         verbose = 2)
tmp <- mbrain; tmp_mat <- multiSVD_obj$laplacian_list$common_laplacian
colnames(tmp_mat) <- paste0("tmp_", 1:ncol(tmp_mat))
set.seed(10); tmp_umap <- Seurat::RunUMAP(tmp_mat)
tmp[["common_laplacian"]] <- Seurat::CreateDimReducObject(tmp_umap@cell.embeddings)
plot1 <- Seurat::DimPlot(tmp, reduction = "common_laplacian",
                        group.by = "label_Savercat", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nCommon Laplacian"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_RNA-geneActivity_common-laplacian.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

#########3

multiSVD_obj <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj,
                                      verbose = 1)
multiSVD_obj <- tiltedCCA:::fine_tuning(input_obj = multiSVD_obj,
                                        verbose = 1)
multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = T,
                                                    bool_modality_2_full = T)

save(multiSVD_obj, mbrain,
     date_of_run, session_info,
     file = "../../../out/main/10x_mouseembryo_tiltedcca_RNA-geneActivity.RData")

########################

set.seed(10)
mbrain[["common_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                        what = "common",
                                                        aligned_umap_assay = "umap",
                                                        seurat_obj = mbrain,
                                                        seurat_assay = "RNA",
                                                        verbose = 1)
set.seed(10)
mbrain[["distinct1_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                           what = "distinct_1",
                                                           aligned_umap_assay = "umap",
                                                           seurat_obj = mbrain,
                                                           seurat_assay = "RNA",
                                                           verbose = 1)
set.seed(10)
mbrain[["distinct2_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                           what = "distinct_2",
                                                           aligned_umap_assay = "umap",
                                                           seurat_obj = mbrain,
                                                           seurat_assay = "RNA",
                                                           verbose = 1)

save(multiSVD_obj, mbrain,
     date_of_run, session_info,
     file = "../../../out/main/10x_mouseembryo_tiltedcca_RNA-geneActivity.RData")

