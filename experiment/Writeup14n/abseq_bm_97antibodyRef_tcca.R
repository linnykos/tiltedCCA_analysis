rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/main/abseq_bm97Ref_preprocessed.RData")
source("bm_97antibodyRef_colorPalette.R")

set.seed(10)

date_of_run <- Sys.time()
session_info <- devtools::session_info()

Seurat::DefaultAssay(bm) <- "RNA"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:30)
bm <- Seurat::FindClusters(bm, resolution = 0.25)

Seurat::DefaultAssay(bm) <- "AB"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:30, reduction = "apca")
bm <- Seurat::FindClusters(bm, resolution = 0.25)

plot1 <-Seurat::DimPlot(bm, reduction = "rna.umap",
                        group.by = "RNA_snn_res.0.25", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nRNA clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14n/abseq_bm97Ref_rna-clustering.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "adt.umap",
                        group.by = "AB_snn_res.0.25", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nADT clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14n/abseq_bm97Ref_adt-clustering.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")

#########

Seurat::DefaultAssay(bm) <- "RNA"
mat_1 <- Matrix::t(bm[["RNA"]]@scale.data[Seurat::VariableFeatures(object = bm),])
Seurat::DefaultAssay(bm) <- "AB"
mat_2 <- Matrix::t(bm[["AB"]]@scale.data)

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
                                            dims_1 = 1:30, dims_2 = 1:30,
                                            center_1 = T, center_2 = T,
                                            normalize_row = T,
                                            normalize_singular_value = T,
                                            recenter_1 = F, recenter_2 = F,
                                            rescale_1 = F, rescale_2 = F,
                                            scale_1 = T, scale_2 = T,
                                            verbose = 1)
multiSVD_obj <- tiltedCCA:::form_metacells(input_obj = multiSVD_obj,
                                           large_clustering_1 = as.factor(bm$RNA_snn_res.0.25), 
                                           large_clustering_2 = as.factor(bm$AB_snn_res.0.25), 
                                           num_metacells = 5000,
                                           verbose = 1)
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 20,
                                         num_neigh = 60,
                                         bool_cosine = T,
                                         bool_intersect = T,
                                         min_deg = 15,
                                         verbose = 2)

tmp <- bm; tmp_mat <- multiSVD_obj$laplacian_list$common_laplacian
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
tmp[["common_laplacian"]] <- Seurat::CreateDimReducObject(tmp_umap_full, key = "commonLapUMAP")
# tmp[["common_laplacian"]] <- Seurat::CreateDimReducObject(tmp_umap, key = "commonLapUMAP")
plot1 <- Seurat::DimPlot(tmp, reduction = "common_laplacian",
                         group.by = "ct", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nCommon Laplacian"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14n/abseq_bm97Ref_common-laplacian.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")

multiSVD_obj <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj,
                                      verbose = 1)
multiSVD_obj <- tiltedCCA:::fine_tuning(input_obj = multiSVD_obj,
                                        verbose = 1)
multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1)

save(multiSVD_obj, bm,
     date_of_run, session_info,
     file = "../../../../out/Writeup14n/abseq_bm97Ref_tcca.RData")

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

save(multiSVD_obj, bm,
     date_of_run, session_info,
     file = "../../../../out/Writeup14n/abseq_bm97Ref_tcca.RData")
