rm(list=ls())
library(Seurat)
library(Signac)
source("reik_colorPalette.R")

load("../../../../out/main/10x_reik_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

Seurat::DefaultAssay(reik) <- "RNA"
mat_1 <- Matrix::t(reik[["RNA"]]@data[Seurat::VariableFeatures(object = reik),])
Seurat::DefaultAssay(reik) <- "ATAC"
mat_2 <- Matrix::t(reik[["ATAC"]]@data[Seurat::VariableFeatures(object = reik),])

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
                                           large_clustering_1 = NULL, 
                                           large_clustering_2 = NULL, 
                                           num_metacells = 5000,
                                           verbose = 1)
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 20,
                                         num_neigh = 100,
                                         bool_cosine = T,
                                         bool_intersect = T,
                                         min_deg = 100,
                                         verbose = 2)
tmp <- reik; tmp_mat <- multiSVD_obj$laplacian_list$common_laplacian
colnames(tmp_mat) <- paste0("tmp_", 1:ncol(tmp_mat))
set.seed(10); tmp_umap <- Seurat::RunUMAP(tmp_mat)@cell.embeddings
tmp_umap_full <- matrix(NA, nrow = ncol(tmp), ncol = 2)
for(i in 1:length(multiSVD_obj$metacell_obj$metacell_clustering_list)){
  idx <- multiSVD_obj$metacell_obj$metacell_clustering_list[[i]]
  tmp_umap_full[idx,] <- rep(tmp_umap[i,], each = length(idx))
}
set.seed(10)
tmp_umap_full <- jitter(tmp_umap_full)
rownames(tmp_umap_full) <- colnames(tmp)
tmp[["common_laplacian"]] <- Seurat::CreateDimReducObject(tmp_umap_full, key = "commonLapUMAP",
                                                          assay = "RNA")
plot1 <- Seurat::DimPlot(tmp, reduction = "common_laplacian",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nCommon Laplacian"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_reik_umap_common-laplacian_v3.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")

#########3

multiSVD_obj <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj,
                                      verbose = 1)
multiSVD_obj <- tiltedCCA:::fine_tuning(input_obj = multiSVD_obj,
                                        verbose = 1)
multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = T,
                                                    bool_modality_2_full = F)

save(multiSVD_obj, reik, peak_granges,
     date_of_run, session_info,
     file = "../../../out/main/10x_reik_tiltedcca.RData")

########################

set.seed(10)
reik[["common_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                      what = "common",
                                                      aligned_umap_assay = "umap",
                                                      seurat_obj = reik,
                                                      seurat_assay = "RNA",
                                                      verbose = 1)
set.seed(10)
reik[["distinct1_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                         what = "distinct_1",
                                                         aligned_umap_assay = "umap",
                                                         seurat_obj = reik,
                                                         seurat_assay = "RNA",
                                                         verbose = 1)
set.seed(10)
reik[["distinct2_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                         what = "distinct_2",
                                                         aligned_umap_assay = "umap",
                                                         seurat_obj = reik,
                                                         seurat_assay = "RNA",
                                                         verbose = 1)

save(multiSVD_obj, reik, peak_granges,
     date_of_run, session_info,
     file = "../../../../out/Writeup14p/Writeup14p_10x_reik_tiltedcca_v3.RData")
